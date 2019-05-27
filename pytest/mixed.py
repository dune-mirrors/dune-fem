import dune.fem
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
import dune.create as create
from dune.ufl import DirichletBC, Space

## simple example of mixed BCs (u[0] with zero Neumann, u[1] with Dirichlet = 2)
dune.fem.parameter.append("parameter")
dune.fem.parameter.append( {"fem.verboserank": -1} )

newtonParameter = {"tolerance": 1e-5, "verbose": "true",
                   "linear.absolutetol": 1e-6, "linear.reductiontol": 1e-6,
                   "linear.preconditioning.method": "ilu",
                   "linear.preconditioning.iterations": 1, "linear.preconditioning.relaxation": 1.2,
                   "linear.verbose": "true"}

dimDomain = 2
dimRange = 2
grid = create.grid("ALUConform", "../data/mixed.dgf", dimgrid=2)

from ufl import TestFunction, TrialFunction, SpatialCoordinate, conditional, dot
uflSpace = Space(dimDomain, dimRange)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

from math import pi,log,sqrt
from ufl import cos, sin,as_vector, dx, ds, grad, inner
exact = as_vector( [sin(3*pi*x[0])*x[1]*x[1], x[0]*x[0]*cos(3.*pi*x[1])] )
laplace = lambda u : grad(grad(u))[0,0]+grad(grad(u))[1,1]
f = as_vector([ -laplace(exact[0])+exact[0], -laplace(exact[1])+exact[1] ])
left = conditional(x[0]<1e-5,1,0)
equation = inner(grad(u), grad(v)) * dx + inner(u,v) * dx - inner(f,v) * dx \
           + 3*pi*x[1]**2*v[0] * left*ds == 0

# spc = create.space("dgonb", grid, dimRange=dimRange, order=1)
spc = create.space("lagrange", grid, dimRange=dimRange, order=1)

def test(operator):
    model = [equation,
            DirichletBC(uflSpace, [None,x[0]**2], 2),       # bottom
            DirichletBC(uflSpace, [exact[0],None], 3),      # top
            DirichletBC(uflSpace, [None,None], 4),          # left
            DirichletBC(uflSpace, exact, 1)]                # right
    parameters={"newton." + k: v for k, v in newtonParameter.items()}

    scheme = create.scheme(operator, model, spc, parameters=parameters)
    solution = spc.interpolate([0,0],name="solution")
    scheme.solve(target=solution)
    l2errA = sqrt( integrate(grid, (solution-exact)**2, 5) )
    grid.hierarchicalGrid.globalRefine(2)
    # note: without the `clear` the code can fail since the new dofs in 'solution' can be nan
    solution.clear()
    scheme.solve(target=solution)
    l2errB = sqrt( integrate(grid, (solution-exact)**2, 5) )
    return solution, l2errA,l2errB

solution, l2errA, l2errB = test("galerkin")
l2eoc = log(l2errA/l2errB)/log(2.)
# print(l2errA,l2errB,l2eoc)
# grid.writeVTK("mixedGalerkin", pointdata={"solution":solution, "exact":exact})
assert abs(l2eoc - (spc.order+1)) < 0.3

solution, l2errA, l2errB = test("h1")
l2eoc = log(l2errA/l2errB)/log(2.)
# print(l2errA,l2errB,l2eoc)
# grid.writeVTK("mixedh1", pointdata={"solution":solution, "exact":exact})
assert abs(l2eoc - (spc.order+1)) < 0.3
