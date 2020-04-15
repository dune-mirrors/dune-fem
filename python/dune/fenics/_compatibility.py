from __future__ import print_function

import dune.common
from dune.ufl import *
from ufl import *
from ufl.log import UFLException
from dune.ufl import DirichletBC, Constant
from dune.ufl import cell as cell_
from dune.fem.function import uflFunction

# return a structured mesh for the unit square
def Point(x,y=None,z=None):
    if z is None:
        if y is None:
            return [x]
        else:
            return [x,y]
    else:
        return [x,y,z]
def BoxMesh( lower, upper, xspacing, yspacing = None, zspaceing = None, **unused ):
    from dune.grid import structuredGrid
    if( zspaceing is None ):
        if( yspacing is None ):
            return structuredGrid(lower, upper, [xspacing] )
        else:
            return structuredGrid(lower, upper, [xspacing,yspacing] )
    else:
        return structuredGrid(lower, upper, [xspacing,yspacing,zspaceing] )
def RectangleMesh( lower, upper, xspacing, yspacing = None, zspaceing = None, **unused ):
    from dune.grid import structuredGrid
    if( zspaceing is None ):
        if( yspacing is None ):
            return structuredGrid(lower, upper, [xspacing] )
        else:
            return structuredGrid(lower, upper, [xspacing,yspacing] )
    else:
        return structuredGrid(lower, upper, [xspacing,yspacing,zspaceing] )
def UnitSquareMesh( xspacing, yspacing = None, zspaceing = None, **unused ):
    if( zspaceing is None ):
        if( yspacing is None ):
            return RectangleMesh([0],[1],xspacing)
        else:
            return RectangleMesh([0,0],[1,1],xspacing,yspacing)
    else:
        return RectangleMesh([0,0,0],[1,1,1],xspacing,yspacing,zspacing)

# create a discrete functions space given a grid view (called mesh in Fenics)
def FunctionSpace( mesh, family, degree=1, dimrange=None, **kwargs ):
    from dune.fem.space import lagrange, dgonb, bdm, raviartThomas
    if( family in ['P','Lagrange', 'CG'] ):
        return lagrange(gridView=mesh, order=degree, dimRange=dimrange, **kwargs )
    elif( family == 'DG' ):
        return dgonb(gridView=mesh, order=degree, dimRange=dimrange, **kwargs )
    elif( family == 'BDM' ):
        return bdm(gridView=mesh, order=degree, dimRange=dimrange, **kwargs )
    elif( family == 'RT' ):
        return raviarThomas(gridView=mesh, order=degree, dimRange=dimrange, **kwargs )
    else:
        raise ValueError('Space with identifier',spacetype,' not known\n')

def VectorFunctionSpace( mesh, family, degree=1, dimrange=None, **kwargs ):
    if dimrange is None:
        dimrange = mesh.dimension
    return FunctionSpace(mesh, family, degree, dimrange, **kwargs)

# creates a discrete function given a discrete space
_counter = -1
def Function( discreteSpace, name='U', **unused ):
    global _counter
    _counter += 1
    return discreteSpace.interpolate([0,]*discreteSpace.dimRange,
            name+"_"+str(_counter))
def interpolate(expr, discreteSpace, name="U", **unuser ):
    global _counter
    _counter += 1
    return discreteSpace.interpolate(expr, name+"_"+str(_counter))

SpatialCoordinate_ = SpatialCoordinate
def SpatialCoordinate(cell,**unused):
    return SpatialCoordinate_(cell_(cell))
    try:
        return SpatialCoordinate_(cell)
    except UFLException:
        return SpatialCoordinate_(cell_(cell))


# solve given equation using galerkin scheme and integrands model
def solve( equation, target, bc = None, solver=None, **kwargs):
    from dune.fem.scheme import galerkin
    from dune.common.checkconfiguration import assertHave, ConfigurationError
    if solver is None:
        try:
            assertHave("HAVE_UMFPACK")
            solver = ("suitesparse","umfpack")
        except ConfigurationError:
            solver = "gmres"

    if bc is None:
        problem = equation
    else:
        try:
            problem = [equation, *bc]
        except TypeError:
            problem = [equation, bc]

    scheme = galerkin( problem, target.space, solver, **kwargs)
    scheme.solve( target=target )

# plot data or grid
def plot( obj, **unused ):
    from dune.fem.plotting import plotPointData as plot
    plot( obj, block=False )
def interactive():
    from matplotlib import pyplot
    pyplot.show()

def errornorm( a, b, normid='L2', **kwargs ):
    from math import sqrt
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
        order += b.order
    except AttributeError:
        pass
    if normid == 'L2':
        error = inner(a - b, a - b)
        return sqrt( integrate(grid,error,2*order+1) )
    else:
        raise ValueError('errornorm with identifier',normid,' not known\n')


def Expression(cpp_code=None, name=None, degree=None, mesh=None, dimRange=None, **kwargs):
    if name is None:
        global _counter
        _counter += 1
        name = "expr"+str(_counter)
    if degree is None:
        degree = kwargs.get("order",1)
    if mesh is None:
        mesh = kwargs.get("cell",None)
    if mesh is None:
        raise ValueError("no mesh provided")
    if type(cpp_code) in [tuple,list]:
        assert dimRange is None or dimRagne==len(cpp_code)
        dimRange = len(cpp_code)
    else:
        cpp_code = [cpp_code]
        dimRange = 0

    uflList = ['acos', 'asin', 'atan', 'cos',
               'cosh', 'exp',
               'pi', 'sin', 'sinh', 'sqrt',
               'tan', 'tanh']
    # creating a dictionary of safe methods
    uflDict = dict([(k, globals()[k]) for k in uflList])
    uflDict["x"] = SpatialCoordinate(mesh)
    for c in kwargs:
        if not c in ["order","cell"]:
            if not isinstance(kwargs[c], Coefficient) and\
               not isinstance(kwargs[c], Indexed):
                if hasattr(ufl.Coefficient,c):
                    raise AttributeError("can not name a constant "+c\
                            +" this will lead to conflicts with existing"\
                            +" attributes on the ufl.Coefficient class")
                uflDict[c] = Constant(kwargs[c],name=c)
            else:
                uflDict[c] = kwargs[c]
    expr = [None,]*max(dimRange,1)
    for i in range(max(dimRange,1)):
        # in ufl M_PI is just pi. This needs to be replaced here
        if cpp_code[i].find('M_PI'):
            cpp_code[i] = cpp_code[i].replace('M_PI', 'pi' )
        expr[i] = eval(cpp_code[i], {}, uflDict)
        # expr[i] = eval(cpp_code[i], {"__builtins__":None}, uflDict)
    # fails with 'x' undefined? expr = [ eval(code) for code in cpp_code ]
    if dimRange == 0:
        func = uflFunction(mesh, name, degree, expr[0], scalar=dimRange==0)
    else:
        func = uflFunction(mesh, name, degree, expr, scalar=dimRange==0)
    # for c in kwargs:
    #     if not c in ["order","cell"]:
    #         if not isinstance(kwargs[c], Coefficient) and\
    #            not isinstance(kwargs[c], Indexed):
    #                assert hasattr(type(func.__impl__),c)
    #                getattr(type(func.__impl__), c).fset(func.__impl__, kwargs[c])
    return func
