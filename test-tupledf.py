from __future__ import print_function, division

from ufl import as_vector, dot, grad, cos, pi, SpatialCoordinate, triangle
from dune.grid import structuredGrid, gridFunction
from dune.fem.space import lagrange, combined, product

x     = SpatialCoordinate(triangle)
exact = as_vector([cos(2.*pi*x[0])*cos(2.*pi*x[1]),dot(x,x)])

grid     = structuredGrid([0, 0], [1, 1], [16, 16])
spc1     = lagrange(grid, dimRange=1, order=1)
spc2     = lagrange(grid, dimRange=1, order=2)
test1    = spc1.interpolate(exact[0],name="test")
test2    = spc2.interpolate(exact[1],name="test")
spc      = combined( spc1, spc2 )
solution = spc.interpolate(exact, name="solution")

space = product( spc1, spc2, components=["p","s"] )
df    = space.interpolate( exact, name="df")
# print(df.dofVector.size,solution.dofVector.size,
#       df.components[0].dofVector.size,df.p.dofVector.size,test1.dofVector.size)
assert df.components[0].dofVector.size == test1.dofVector.size
assert df.s.dofVector.size == test2.dofVector.size
assert df.dofVector.size   == solution.dofVector.size
df.interpolate(solution)
solution.interpolate(df)
test1.interpolate(df.p)
df.s.interpolate(test2)
df.components[0].interpolate(solution[0])
df.p.interpolate(solution[0])
