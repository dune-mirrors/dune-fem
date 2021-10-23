from mpi4py import MPI

import numpy, math, sys
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.grid import structuredGrid as leafGridView
from dune.fem import threading, parameter
from dune.fem.space import dgonb as dgSpace # dglegendre as dgSpace
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.scheme import molGalerkin as molSolutionScheme
from dune.fem.function import integrate, uflFunction
from dune.ufl import Constant, DirichletBC
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal
from ufl import dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional
from ufl import as_vector, avg, jump, dS, CellVolume, FacetArea, atan, tanh, sin

parameter.append({"fem.verboserank": 0})

def compute(space,epsilon,weakBnd,skeleton, mol=None):
    u    = TrialFunction(space)
    v    = TestFunction(space)
    n    = FacetNormal(space)
    he   = avg( CellVolume(space) ) / FacetArea(space)
    hbnd = CellVolume(space) / FacetArea(space)
    x    = SpatialCoordinate(space)

    exact = uflFunction( space.gridView, name="exact", order=3, ufl=sin(x[0]*x[1]))
    uh = space.interpolate(exact,name="solution")

    # diffusion factor
    eps = Constant(epsilon,"eps")
    # transport direction and upwind flux
    b    = as_vector([1,0])
    hatb = (dot(b, n) + abs(dot(b, n)))/2.0
    # characteristic function for left/right boundary
    dD   = conditional((1+x[0])*(1-x[0])<1e-10,1,0)
    # penalty parameter
    beta = Constant( 20*space.order**2,"beta")

    rhs           = -( div(eps*grad(exact)-b*exact) ) * v  * dx
    aInternal     = dot(eps*grad(u) - b*u, grad(v)) * dx
    aInternal    -= eps*dot(grad(exact),n)*v*(1-dD)*ds

    diffSkeleton  = eps*beta/he*jump(u)*jump(v)*dS -\
                    eps*dot(avg(grad(u)),n('+'))*jump(v)*dS -\
                    eps*jump(u)*dot(avg(grad(v)),n('+'))*dS
    if weakBnd:
        diffSkeleton += eps*beta/hbnd*(u-exact)*v*dD*ds -\
                        eps*dot(grad(exact),n)*v*dD*ds
    advSkeleton   = jump(hatb*u)*jump(v)*dS
    if weakBnd:
        advSkeleton  += ( hatb*u + (dot(b,n)-hatb)*exact )*v*dD*ds

    if skeleton:
        form          = aInternal + diffSkeleton + advSkeleton
    else:
        form          = aInternal

    if weakBnd and skeleton:
        strongBC = None
    else:
        strongBC = DirichletBC(space,exact,dD)

    if space.storage[0] == "numpy":
        solver={"solver":("suitesparse","umfpack"),
                "parameters":{"newton.verbose": True,
                              "newton.linear.verbose": False,
                              "newton.linear.tolerance":1e-5,
                }
               }
    else:
        solver={"solver":"bicgstab",
                "parameters":{"newton.linear.preconditioning.method":"ilu",
                              "newton.linear.tolerance":1e-13,
                              "newton.verbose": True, "newton.linear.verbose": False}
               }
    if mol == 'mol':
        scheme = molSolutionScheme([form==rhs,strongBC], **solver)
    else:
        scheme = solutionScheme([form==rhs,strongBC], **solver)

    eoc = []
    info = scheme.solve(target=uh)

    error = dot(uh-exact,uh-exact)
    error0 = math.sqrt( integrate(gridView,error,order=5) )
    print(error0," # output",flush=True)
    for i in range(3):
        gridView.hierarchicalGrid.globalRefine(1)
        uh.interpolate(exact)
        scheme.solve(target=uh)
        error = dot(uh-exact,uh-exact)
        error1 = math.sqrt( integrate(gridView,error,order=5) )
        eoc   += [ math.log(error1/error0) / math.log(0.5) ]
        print(i,error0,error1,eoc," # output",flush=True)
        error0 = error1

    # print(space.order,epsilon,eoc)
    if (eoc[-1]-(space.order+1)) < -0.1:
        print("ERROR:",space.order,epsilon,eoc)
    return eoc


storage = "istl"
threading.use = 4

def newGridView():
    return leafGridView([-1, -1], [1, 1], [4, 4])

test = True
for i in range(10):
  print(i,"dgSpace, 1e-5, True, True")
  gridView = newGridView()
  space    = dgSpace(gridView, order=2, storage=storage)
  eoc = compute(space,1e-5,True,True) # , 'mol')
  test = test and (eoc[-1]-(space.order+1)) > -0.1

  print(i,"dgSpace, 1, True, True")
  gridView = newGridView()
  space    = dgSpace(gridView, order=2, storage=storage)
  eoc = compute(space,1,True,True)
  test = test and (eoc[-1]-(space.order+1)) > -0.1

  print(i,"lagrange, 1, True, True")
  gridView = newGridView()
  space    = lagrange(gridView, order=2, storage=storage)
  eoc = compute(space,1,True,True)
  test = test and (eoc[-1]-(space.order+1)) > -0.1

  print(i,"lagrange, 1, False, True")
  gridView = newGridView()
  space    = lagrange(gridView, order=2, storage=storage)
  eoc = compute(space,1,False,True)
  test = test and (eoc[-1]-(space.order+1)) > -0.1

  print(i,"lagrange, 1, False, False")
  gridView = newGridView()
  space    = lagrange(gridView, order=2, storage=storage)
  eoc = compute(space,1,False,False)
  test = test and (eoc[-1]-(space.order+1)) > -0.1

  if not test: break

print("-----------------------")
sys.exit(0 if test else 1)
