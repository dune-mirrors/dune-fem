from __future__ import print_function

import dune.common
from dune.ufl import *
from ufl import *
from ufl.log import UFLException
from dune.ufl import DirichletBC, NamedConstant
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

# convenience function for Fenics' Constant
def Constant( value, **unused ):
    try:
        if len(value)>0:
            return as_vector(value)
    except:
        pass
    return value

# create a discrete functions space given a grid view (called mesh in Fenics)
def FunctionSpace( mesh, family, degree=1, dimrange=None, **kwargs ):
    from dune.fem.space import lagrange, dgonb, bdm, raviartThomas
    if( family in ['P','Lagrange', 'CG'] ):
        return lagrange(view=mesh, order=degree, dimrange=dimrange, **kwargs )
    elif( family == 'DG' ):
        return dgonb(view=mesh, order=degree, dimrange=dimrange, **kwargs )
    elif( family == 'BDM' ):
        return bdm(view=mesh, order=degree, dimrange=dimrange, **kwargs )
    elif( family == 'RT' ):
        return raviarThomas(view=mesh, order=degree, dimrange=dimrange, **kwargs )
    else:
        raise ValueError('Space with identifier',spacetype,' not known\n')

def VectorFunctionSpace( mesh, family, degree=1, dimrange=None, **kwargs ):
    if dimRange is None:
        dimRange = mesh.dimension
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
def solve( equation, target, bc = None, solver='gmres', **kwargs):
    from dune.fem.scheme import galerkin
    if bc is None:
        problem = equation
    else:
        try:
            problem = [equation, *bc]
        except TypeError:
            problem = [equation, bc]
    scheme = galerkin( problem, target.space, solver, **kwargs)
            # parameters={"newton.verbose":True,"newton.linear.verbose":True})
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
        return sqrt( integrate(grid,error,2*order+1)[0] )
    else:
        raise ValueError('errornorm with identifier',normid,' not known\n')

def Expression(cpp_code=None, name="tmp", order=None, view=None, **kwargs):
    if order is None:
        order = kwargs.get("degree",1)
    if view is None:
        view = kwargs.get("cell",None)
    if view is None:
        raise ValueError("no view provided")
    if type(cpp_code) in [tuple,list]:
        dimRange = len(cpp_code)
    else:
        cpp_code = [cpp_code]
        dimRange = 0

    x = SpatialCoordinate(view)
    coeffs = {}
    namedCoeffs = []
    for c in kwargs:
        if not c in ["degree","cell"]:
            namedCoeffs.append(NamedConstant(view,c))
            coeffs[c] = "namedCoeffs["+str(len(namedCoeffs)-1)+"]"
    if len(coeffs)>0:
        import re
        coeffs = dict((re.escape(k), v) for k, v in coeffs.items())
        pattern = re.compile("|".join(coeffs.keys()))
        cpp_code = [ pattern.sub(lambda m: coeffs[re.escape(m.group(0))], code) for code in cpp_code ]
    expr = [None,]*max(dimRange,1)
    for i in range(max(dimRange,1)):
        expr[i] = eval(cpp_code[i])
    # fails with 'x' undefined? expr = [ eval(code) for code in cpp_code ]
    func = uflFunction(view, name, order, expr)
    for c in namedCoeffs:
        func.__dict__[c.name] = kwargs[c.name]
    return func[0] if dimRange==0 else func
