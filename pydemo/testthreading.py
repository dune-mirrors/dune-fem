from mpi4py import MPI

import numpy, math, time
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.grid import structuredGrid as leafGridView
# from dune.alugrid import aluCubeGrid as leafGridView
import dune.fem
from dune.grid import cartesianDomain
from dune.fem.space import dgonb as dgSpace # dglegendre as dgSpace
from dune.fem.space import lagrange
# from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.scheme import molGalerkin as solutionScheme
from dune.fem.function import integrate, uflFunction
from dune.fem.operator import linear
from dune.ufl import Constant, DirichletBC
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal
from ufl import dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional
from ufl import as_vector, avg, jump, dS, CellVolume, FacetArea, atan, tanh, sin

def model(space,epsilon,weakBnd,skeleton):
    u    = TrialFunction(space)
    v    = TestFunction(space)
    n    = FacetNormal(space)
    he   = avg( CellVolume(space) ) / FacetArea(space)
    hbnd = CellVolume(space) / FacetArea(space)
    x    = SpatialCoordinate(space)

    #exact = sin(x[0]*x[1]) # atan(1*x[1])
    exact = uflFunction( space.grid, name="exact", order=3, ufl=sin(x[0]*x[1]))

    # diffusion factor
    eps = 1 # Constant(epsilon,"eps")
    # transport direction and upwind flux
    b    = as_vector([1,0])
    hatb = (dot(b, n) + abs(dot(b, n)))/2.0
    # characteristic function for left/right boundary
    dD   = conditional((1+x[0])*(1-x[0])<1e-10,1,0)
    # penalty parameter
    beta = Constant( 10*space.order**2 if space.order > 0 else 1,"beta")

    rhs           = ( -div(eps*grad(exact)-b*exact) + exact) * v  * dx
    aInternal     = (dot(eps*grad(u) - b*u, grad(v)) + dot(u,v)) * dx
    diffSkeleton  = eps*beta/he*jump(u)*jump(v)*dS -\
                    eps*dot(avg(grad(u)),n('+'))*jump(v)*dS -\
                    eps*jump(u)*dot(avg(grad(v)),n('+'))*dS
    diffSkeleton -= eps*dot(grad(exact),n)*v*(1-dD)*ds
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
        strongBC = None # DirichletBC(space,exact,dD)

    if space.storage[0] == "fem":
        solver={"solver":("suitesparse","umfpack")}
    else:
        solver={"solver":"bicgstab",
                "parameters":{"newton.linear.preconditioning.method":"ilu",
                              "newton.linear.tolerance":1e-13}
               }
    scheme = solutionScheme([form==rhs,strongBC], **solver)
    uh = space.interpolate([0],name="solution")
    A = linear(scheme)
    return scheme, uh, A

def compute(scheme, uh, A ):
    start = time.time()
    scheme(uh.copy(),uh)
    runTime = [time.time()-start]
    start = time.time()
    scheme.jacobian(uh,A)
    runTime += [time.time()-start]
    return runTime

storage = "fem"

def newGridView(N=4):
    return leafGridView( [-1, -1], [1, 1], [N, N ])
    ctor = cartesianDomain( [-1, -1], [1, 1], [N, N ])
    return leafGridView(ctor)

def test(spaceCtor,skeleton):
    defaultThreads = dune.fem.threading.use
    # the operator is run once when setting up the linear operator in the 'model'
    gridView = newGridView(4)
    space    = spaceCtor(gridView, order=2, storage=storage)
    scheme, uh, A = model(space,1,True,skeleton)
    compute(scheme,uh,A)

    gridView = newGridView(N=400)
    space    = spaceCtor(gridView, order=2, storage=storage)
    scheme, uh, A = model(space,1,True,skeleton)

    print("---------------------")

    # time with the default number of threads (1 if no environment variable is set)
    runTime = compute(scheme,uh,A)
    print(dune.fem.threading.use," thread used: ",runTime,flush=True)

    # time with 2 threads
    dune.fem.threading.use = 2
    runTime = compute(scheme,uh,A)
    print(dune.fem.threading.use," threads used: ",runTime,flush=True)

    # time with 4 threads
    dune.fem.threading.use = 4
    runTime = compute(scheme,uh,A)
    print(dune.fem.threading.use," threads used: ",runTime,flush=True)

    # time with max number of threads
    dune.fem.threading.use = 8
    runTime = compute(scheme,uh,A)
    print(dune.fem.threading.use," threads used: ",runTime,flush=True)

    # time with max number of threads
    dune.fem.threading.use = 16
    runTime = compute(scheme,uh,A)
    print(dune.fem.threading.use," threads used: ",runTime,flush=True)

    dune.fem.threading.use = defaultThreads

print("DGSpace:")
test(dgSpace,True)
print("Lagrange (with skeleton):")
test(lagrange,True)
print("Lagrange (no skeleton):")
test(lagrange,False)
