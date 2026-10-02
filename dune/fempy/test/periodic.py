# Periodic test: we solve grad(u).grad(v) + uv = fv
# with f so that u(x,y)=sin(x pi/4) is the exact solution on Omega=(0,16)x(0,4).
#
# The solution has 0 Neumann b.c. at the top and bottom but not on the left and right
# so without periodicity in x addition flux b.c. would be needed.
# As a result the EOC is clearly wrong if no periodicity is used.
#
# Tested with a ALUConform and Lagrange space with order 1 and 4.
#
# Note: using u(x,y)=cos(x pi/4) works also without periodicity so can be used as a test.
#
# The given dgfStr should produce the right periodic ALU grid.
# Using dune.periodic.flatPeriodicGrid does work as expected - after some fixes are applied
# to dune-grid and dune-alugrid (branch ...). The fix in dune-alugrid is not perfect
# requiring that alugrid_assert( hsfc.find( sfc.index( center ) ) == hsfc.end() );
# is removed in dune/alugrid/3d/gridfactory.cc.
#
# An extended test a dune.fem.globalRefine can be used with prolongation.
# The error is compute before and after a globalrefine(1,[uh])
# and we check that the error is the same.
# Note: for some reason the dune.fem.globalRefine does not lead to anything happening with dune.periodic.

dgfStr = """\
DGF
Interval
0.0 0.0
16.0 4.0
90 30
#
BOUNDARYDOMAIN
3 0.0 0.0 16.0 0.0
4 0.0 4.0 16.0 4.0
default 1
#
"""
dgfPer = """
PERIODICFACETRANSFORMATION
1 0, 0 1 + 16 0
#
"""
dgfStr += dgfPer

import math, sys, io
import numpy as np
from dune.generator import algorithm
from dune.grid import reader
from dune.alugrid import aluConformGrid as hostGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
from dune.fem.function import gridFunction
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem import integrate, globalRefine
from dune.ufl import Constant
from ufl import ( TestFunction, TrialFunction, SpatialCoordinate, FacetNormal,
                  pi, dx, ds, grad, div, grad, dot, inner, sqrt, exp, atan, cos, sin, conditional,
                  as_vector, avg, jump, dS, CellVolume, FacetArea )


def run(space, exact):
    gridView = space.gridView
    print(f"vertices={gridView.size(2)} and dofs={space.size}")

    u,v,x  = TrialFunction(space), TestFunction(space), SpatialCoordinate(space)
    f = -div(grad(exact)) + exact
    g    = exact
    aInternal    = ( dot(grad(u), grad(v)) + u*v - f*v ) * dx
    form         = aInternal
    scheme = solutionScheme(form==0, solver="cg",
                            parameters={"linear.preconditioning.method":"ssor"}
                           )
    uh = space.function(name="solution")
    info = scheme.solve(target=uh)
    return uh

domain = (reader.dgfString, dgfStr)
gridView = hostGridView(domain, dimgrid=2)

# this approach works - but requires some fixes:
# 1. in dune-grid to avoid dgfGF reader instantiation
# 2. issue in SFC in dune-alugrid, setting sfc=None in the hostgrid does not
#             help since the newly created grid in dune-periodic is the issue.
# from dune.periodic import flatPeriodicGrid as pGrid
# gridView = pGrid(gridView, (1,0))

gridView = adaptiveGridView( gridView )
space = solutionSpace(gridView, order=1)
x = SpatialCoordinate(space)
# switching to cos here works fine due to Neumann zero boundary value
# exact = cos(x[0]*pi/4)
exact = sin(x[0]*pi/4)

levels = 2
for i in range(levels):
    uh = run(space, exact)
    l2error = dot(uh-exact,uh-exact)
    h1error = dot(grad(uh-exact),grad(uh-exact))
    errors = [np.sqrt(e) for e in integrate([l2error,h1error])]
    print('\t | u_h - u | =', '{:0.5e}'.format(errors[0]))
    print('\t | grad(uh - u) | =', '{:0.5e}'.format(errors[1]))
    if i < levels-1:
        if False: # use dune-fem's global refine and check prolongation
            print("Before grid:",gridView.size(2),gridView.hierarchicalGrid.leafView.size(2))
            globalRefine(2,[uh])
            print("After grid:",gridView.size(2),gridView.hierarchicalGrid.leafView.size(2))
            errors2 = [np.sqrt(e) for e in integrate([l2error,h1error])]
            print('\t | u_h - u | =', '{:0.5e}'.format(errors[0]))
            print('\t | grad(uh - u) | =', '{:0.5e}'.format(errors[1]))
            assert np.all(np.isclose(errors2,errors))
        else:
            gridView.hierarchicalGrid.globalRefine(1)
