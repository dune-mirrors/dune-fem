from __future__ import print_function

import dune.common
from dune.ufl import *
from ufl import *
from ufl.log import UFLException
from dune.ufl import DirichletBC
from dune.ufl import cell as cell_

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

def Constant( value, **unused ):
    try:
        if len(value)>0:
            return as_vector(value)
    except:
        pass
    return value

# create a discrete functions space given a grid view (called mesh in Fenics)
def FunctionSpace( gridview, spacetype, order=1, dimRange=None, **kwargs ):
    from dune.fem.space import lagrange, dgonb
    if( spacetype in ['P','lagrange'] ):
        return lagrange( gridview,order,dimRange,**kwargs )
    elif( spacetype == 'DG' ):
        return dgonb( gridview, order,dimRange,**kwargs )
    else:
        raise ValueError('Space with identifier',spacetype,' not known\n')
def VectorFunctionSpace( gridview, spacetype, order=1, dimRange=None, **kwargs ):
    if dimRange is None:
        dimRange = gridview.dimension
    return FunctionSpace(gridview, spacetype, order, dimRange, **kwargs)

# creates a discrete function given a discrete space
def Function( discreteSpace, name='U', **unused ):
    return discreteSpace.interpolate([0,]*discreteSpace.dimRange, name)

SpatialCoordinate_ = SpatialCoordinate
def SpatialCoordinate(cell,**unused):
    return SpatialCoordinate_(cell_(cell))
    try:
        return SpatialCoordinate_(cell)
    except UFLException:
        return SpatialCoordinate_(cell_(cell))


# solve given equation using galerkin scheme and integrands model
def solve( equation, solution, bc = None, solver='gmres', **kwargs):
    from dune.fem.scheme import galerkin
    if bc is None:
        problem = equation
    else:
        try:
            problem = [equation, *bc]
        except TypeError:
            problem = [equation, bc]
    scheme = galerkin( problem, solution.space, solver, **kwargs,
            parameters={"newton.verbose":True,"newton.linear.verbose":True})
    scheme.solve( target=solution )

# plot data or grid
def plot( obj, **unused ):
    from dune.fem.plotting import plotPointData as plot
    plot( obj, block=False )
def interactive():
    from matplotlib import pyplot
    pyplot.show()

def errornorm( a, b, normid='L2', **kwargs ):
    from dune.fem.function import integrate
    grid = kwargs.get("grid",None)
    order = kwargs.get("order",0)
    if grid is None:
        try:
            grid = a.grid
        except AttributeError:
            pass
    if grid is None:
        try:
            grid = b.grid
        except AttributeError:
            pass
    if grid is None:
        raise ValueError("can not extract grid from arguments")
    try:
        order += a.order
    except AttributeError:
        pass
    try:
        order += a.order
    except AttributeError:
        pass
    print("norm order=",order)
    if normid == 'L2':
        error = inner(a - b, a - b)
        return integrate(grid,error,2*order)
    else:
        raise ValueError('errornorm with identifier',normid,' not known\n')
