#!/usr/bin/env python3

# <markdowncell>
# # Advection-Diffusion: Discontinuous Galerkin Method with Upwinding
# So far we have been using Lagrange spaces of different order to solve our
# PDE. In the following we show how to use Discontinuous Galerkin method to
# solve an advection dominated advection-diffusion probllem:
# \begin{align*}
# -\varepsilon\triangle u + b\cdot\nabla u &= f
# \end{align*}
# with Dirichlet boundary conditions. Here $\varepsilon$ is a small
# constant and $b$ a given vector.
# <codecell>

from mpi4py import MPI

import numpy, math
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.grid import structuredGrid as leafGridView
from dune.fem.space import dgonb as dgSpace # dglegendre as dgSpace
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.function import integrate
from dune.ufl import Constant, DirichletBC
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal
from ufl import dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional
from ufl import as_vector, avg, jump, dS, CellVolume, FacetArea, atan, tanh, sin

def compute(space,epsilon,weakBnd):
    u    = TrialFunction(space)
    v    = TestFunction(space)
    n    = FacetNormal(space)
    he   = avg( CellVolume(space) ) / FacetArea(space)
    hbnd = CellVolume(space) / FacetArea(space)
    x    = SpatialCoordinate(space)

    exact = sin(x[0]*x[1]) # atan(1*x[1])

    # diffusion factor
    eps = Constant(epsilon,"eps")
    # transport direction and upwind flux
    b    = as_vector([1,0])
    hatb = (dot(b, n) + abs(dot(b, n)))/2.0
    # characteristic function for left/right boundary
    dD   = conditional((1+x[0])*(1-x[0])<1e-10,1,0)
    # penalty parameter
    beta = Constant( 10*space.order**2 if space.order > 0 else 1,"beta")

    rhs           = -( div(eps*grad(exact)-b*exact) ) * v  * dx
    aInternal     = dot(eps*grad(u) - b*u, grad(v)) * dx
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
    form          = aInternal + diffSkeleton + advSkeleton

    if weakBnd:
        strongBC = None
    else:
        strongBC = DirichletBC(space,exact,dD)

    if space.storage[0] == "fem":
        solver={"solver":("suitesparse","umfpack")}
    else:
        solver={"solver":"bicgstab",
                "parameters":{"newton.linear.preconditioning.method":"ilu",
                              "newton.linear.tolerance":1e-13}
               }
    scheme = solutionScheme([form==rhs,strongBC], **solver)

    # <markdowncell>
    # <codecell>

    uh = space.interpolate([0],name="solution")
    eoc = []
    info = scheme.solve(target=uh)
    error0 = math.sqrt( integrate(gridView,dot(uh-exact,uh-exact),order=5) )
    for i in range(3):
        uh.interpolate(0)
        scheme.solve(target=uh)
        error1 = math.sqrt( integrate(gridView,dot(uh-exact,uh-exact),order=5) )
        eoc   += [ math.log(error1/error0) / math.log(0.5) ]
        print(i,error0,error1,eoc)
        error0 = error1
        gridView.hierarchicalGrid.globalRefine(1)

    # print(space.order,epsilon,eoc)
    if (eoc[-1]-(space.order+1)) < -0.1:
        print("ERROR:",space.order,epsilon,eoc)
    assert (eoc[-1]-(space.order+1)) > -0.1
    return eoc

storage = "fem"

def newGridView():
    return leafGridView([-1, -1], [1, 1], [4, 4])

gridView = newGridView()
space    = dgSpace(gridView, order=2, storage=storage)
eoc = compute(space,1e-5,True)

gridView = newGridView()
space    = dgSpace(gridView, order=2, storage=storage)
eoc = compute(space,1,True)

#gridView = newGridView()
#space    = dgSpace(gridView, order=3, storage=storage)
#eoc = compute(space,1e-5,True)

gridView = newGridView()
space    = lagrange(gridView, order=2, storage=storage)
eoc = compute(space,1,True)

gridView = newGridView()
space    = lagrange(gridView, order=2, storage=storage)
eoc = compute(space,1,False)
