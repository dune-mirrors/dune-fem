from __future__ import print_function

import dune.common
from dune.ufl import *
from ufl import *

# return a structured mesh for the unit square
def UnitSquareMesh( xspacing, yspacing = None, zspaceing = None, **unused ):
    from dune.grid import structuredGrid
    if( zspaceing is None ):
        if( yspacing is None ):
            return structuredGrid([0], [1], [xspacing] )
        else:
            return structuredGrid([0,0], [1,1], [xspacing,yspacing] )
    else:
        return structuredGrid([0,0,0], [1,1,1], [xspacing,yspacing,zspaceing] )


# create a discrete functions space given a grid view (called mesh in Fenics)
def FunctionSpace( gridview, spacetype, order=1, dimRange=1, **unused ):
    from dune.fem.space import lagrange, dgonb
    if( spacetype == 'P' ):
        return lagrange( gridview, dimRange, order )
    elif( spacetype == 'DG' ):
        return dgonb( gridview, dimRange, order )
    else:
        print('Space with identifier',spacetype,' not known\n')

# creates a discrete function given a discrete space
def Function( discreteSpace, name='U', **unused ):
    return discreteSpace.interpolate([0], name)

# solve given equation using galerkin scheme and integrands model
def solve( equation, solution, bc = None, solver='cg', **unused):
    from dune.fem.scheme import galerkin
    scheme = galerkin( equation, solution.space, solver )
    scheme.solve( target=solution )

# plot data or grid
def plot( obj, **unused ):
    from dune.fem.plotting import plotPointData as plot
    plot( obj )


def errornorm( a, b, normid='L2',**unused ):
    from dune.fem.function import integrate
    error = a - b
    try:
        return sqrt( error.integrate() )
    except AttributeError:
        print('errornorm failure\n')
