from ufl import *

import math
import dune.fem
from dune.fem.function import integrate
import dune.create as create
from dune.ufl import DirichletBC, Space

dimRange = 12  # needs to be >= 4, test with 4,8,11
grid = create.grid("ALUConform", "../data/mixed.dgf", dimgrid=2)

from ufl import SpatialCoordinate
uflSpace = dune.ufl.Space(2, dimRange)
x = SpatialCoordinate(uflSpace.cell())
from math import pi,log,sqrt
from ufl import cos,sin,as_vector
exact = as_vector( [sin(3*pi*x[0]), x[1]*x[1], x[0]*x[0], cos(3.*pi*x[1])]+[0]*(dimRange-4) )
v1 = integrate(grid, exact, 5).two_norm

space = create.space("Lagrange", grid, dimRange=dimRange, order=1)
u = space.interpolate(exact,name="u")
v2 = integrate(grid, u, 5).two_norm

print(v1,v2,v1-v2)
v3 = integrate(grid, inner(grad(u[11]),grad(u[11])), 5)
