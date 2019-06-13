import time
from dune.grid import structuredGrid
from dune.fem import parameter
from dune.fem.function import integrate
import dune.create as create
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, exp,\
                dx, grad, inner, as_vector, replace, sqrt, dot,\
                derivative, action
from ufl.algorithms.apply_derivatives import apply_derivatives
from dune.fem.operator import linear as linearOperator

parameter.append({"fem.verboserank": -1})

grid = structuredGrid([0, 0], [1, 1], [40, 40])

order     = 2
dimR      = 2
spaceName = "lagrange"
space = create.space(spaceName, grid, dimRange=dimR, order=order)

arg   = space.interpolate(as_vector([1,]*dimR), name='arg')
destA = space.interpolate([0,]*dimR, name='destA')              # note: this uses an efficient approach in C++
destB = space.interpolate(as_vector([0.,]*dimR), name='destB')   # note: this needs to generate a ufl local function in the current implementation
destC = space.interpolate(as_vector([0,]*dimR), name='destC')
destD = space.interpolate(as_vector([0,]*dimR), name='destD')
destE = space.interpolate(as_vector([0,]*dimR), name='destE')

u  = TrialFunction(space)
v  = TestFunction(space)
x = SpatialCoordinate(space.cell())

ubar = space.interpolate(as_vector([dot(x,x),]*dimR),name="ubar")
a    = ( inner(0.5*dot(u,u), v[0]) +
         inner(u[0]*grad(u), grad(v)) ) * dx
op   = create.operator("galerkin", a, space)
A    = linearOperator(op)
op.jacobian(ubar,A)
A(arg,destA)

da   = apply_derivatives(derivative(action(a, ubar), ubar, u))
dop  = create.operator("galerkin", da, space)
dop(arg,destB)
err = integrate(grid, (destA-destB)**2, 5)
# print("error=",err)
assert(err < 1e-15)

A = linearOperator(dop)
dop.jacobian(arg,A)
A(arg,destC)
err = integrate(grid, (destA-destC)**2, 5)
# print("error=",err)
assert(err < 1e-15)

###############################################################

op(ubar,destA)

lina  = ( inner(0.5*dot(ubar,u), v[0]) +
          inner(ubar[0]*grad(u), grad(v)) ) * dx
linop = create.operator("galerkin", lina, space)
linop(ubar,destD)
err = integrate(grid, (destA-destD)**2, 5)
# print("error=",err)
assert(err < 1e-15)

A = linearOperator(linop)
linop.jacobian(arg,A)
A(ubar,destE)
err = integrate(grid, (destA-destE)**2, 5)
# print("error=",err)
assert(err < 1e-15)

# grid.writeVTK('optest', pointdata=[destA,destB,destC,destD,destE])
