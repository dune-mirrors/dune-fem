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

def Constant( value, dimRange=1, **unused ):
    return as_vector( [value,] * dimRange )

# create a discrete functions space given a grid view (called mesh in Fenics)
def FunctionSpace( gridview, spacetype, order=1, dimRange=1, **unused ):
    from dune.fem.space import lagrange, dgonb
    if( spacetype == 'P' or spacetype == 'Lagrange' or spacetype == 'CG' ):
        return lagrange( gridview, dimRange, order )
    elif( spacetype == 'DG' ):
        return dgonb( gridview, dimRange, order )
    else:
        print('Space with identifier',spacetype,' not known\n')

# creates a discrete function given a discrete space
def Function( discreteSpace, name='U', **unused ):
    return discreteSpace.interpolate([0], name)

# solve given equation using galerkin scheme and integrands model
def solve( equation, target, bc = None, solver='gmres', **unused):
    from dune.fem.scheme import galerkin
    if( bc is None ):
        scheme = galerkin( equation, target.space, solver )
    else:
        scheme = galerkin( [equation,bc], target.space, solver )
    scheme.solve( target=target )

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
