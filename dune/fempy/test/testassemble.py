from ufl import *
from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem import assemble
from dune.fem.function import gridFunction
g = structuredGrid([0,0],[1,1],[10,10])
s = lagrange(g)
x = SpatialCoordinate(s)

# cam add a grid attribute through dx?
# Or provide a dune.fem.ufl.SpatialCoordinate that takes a gridView and stores it
#    or even a space to have an 'order'.
I = assemble( x[0]*x[1]*dx, gridView=g, order=2 )
print("integral:", I)
# the following works:
f = gridFunction(x[0], gridView=g, order=2, name="x")
I = assemble( f*x[1]*dx )
print("integral:", I)

v = TestFunction(s)
b = assemble( x[0]*x[1]*v*dx )
print("rhs:", b.scalarProductDofs(b))

u = TrialFunction(s)
A = assemble( x[0]*x[1]*u*v*dx )
print("matrix",A.as_numpy.count_nonzero())

from scipy.sparse.linalg import spsolve as solve
A,b = assemble( x[0]*x[1]*u*v*dx == x[0]*x[1]*v*dx )
print(b.as_numpy)
x = solve(A.as_numpy,b.as_numpy[:])
print(x)
