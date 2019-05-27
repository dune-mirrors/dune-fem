from __future__ import print_function, division

from ufl import as_vector, grad, cos, pi, SpatialCoordinate, triangle
from dune.grid import structuredGrid, cartesianDomain
from dune.fem.space import lagrange
from dune.alugrid import aluConformGrid as gridManager

grid = gridManager(cartesianDomain([0,0],[1,1],[16,16]))
# grid  = structuredGrid([0, 0], [1, 1], [16, 16])
space = lagrange(grid, dimRange=1, order=1)
x = SpatialCoordinate(triangle)
exact = as_vector([cos(2.*pi*x[0])*cos(2.*pi*x[1])])
solution = space.interpolate(exact, name="solution")

test = space.interpolate(solution[0], name="tmp")
test = space.interpolate(grad(solution)[0,0], name="tmp")
