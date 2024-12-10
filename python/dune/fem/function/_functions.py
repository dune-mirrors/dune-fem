from __future__ import absolute_import, division, print_function, unicode_literals

import io, types, inspect, logging, functools
logger = logging.getLogger(__name__)
import numpy as np
from numpy import dtype as np_dtype

import dune.ufl
import dune.grid
import dune.fem.space
import dune.models.localfunction

from dune.fem.deprecated import deprecated

def _integrate(grid,expression,order=None):
    try:
        return expression.integrate()
    except AttributeError:
        return _uflFunction(grid,"tmp",order,expression).integrate()
def integrate(grid,expression,order=None):
    deprecated("dune.fem.function.integrate is deprecated use dune.fem.integrate instead. New signature is (expr, gridView, order)")
    return _integrate(grid,expression,order)

# used to convert python functions with (x) arguments into gridfunctions -
# this is not yet handles by the 'gridFunction' function
def _globalFunction(gridView, name, order, value):
    gf = gridView.function(value,name=name,order=order)
    return dune.ufl.GridFunction(gf)
def globalFunction(gridView, name, order, value):
    deprecated("dune.fem.function.globalFunction is deprecated use dune.fem.function.gridFunction instead. New signature is (expr, gridView, order,order)")
    return _globalFunction(gridView, name, order, value)

# used to convert python functions with (e,x) arguments into gridfunctions -
# this is not yet handles by the 'gridFunction' function
def _localFunction(gridView, name, order, value):
    gf = gridView.function(value,name=name,order=order)
    return dune.ufl.GridFunction(gf)
def localFunction(gridView, name, order, value):
    deprecated("dune.fem.function.localFunction is deprecated use dune.fem.function.gridFunction instead. New signature is (expr, gridView, order,order)")
    return _localFunction(gridView, name, order, value)

def gridFunction(expr=None,gridView=None,*,name=None,order=None, fctName=None, view=None, **kwargs):
    from ufl.core.expr import Expr
    from ufl import as_vector
    if view is not None:
        assert gridView is None
        deprecated("'view' argument changed to 'gridView' in dune.fem.function.gridFunction")
        gridView = view
    if name is None:
        name = "_tmp"

    # first check if 'expr' is already a grid function wrapper
    try:
        if expr.gridView is not None:
            return expr
    except AttributeError:
        pass
    # next change a list/tuple into a ufl expression
    if isinstance(expr, list) or isinstance(expr, tuple):
        expr = as_vector(expr)
    # change a constant into a ufl expression
    if isinstance(expr,int) or isinstance(expr,float):
        expr = dune.ufl.Constant(expr)

    # if 'expr' is a ufl expression use the uflFunction
    if isinstance(expr, Expr): # use the old uflFunction
        return _uflFunction(gridView=gridView,name=name,order=order,ufl=expr)
    # this could be a python functon
    if isinstance(expr, types.FunctionType):
        gl = len(inspect.getfullargspec(expr)[0])
        if gl == 1:   # global function
            func = _globalFunction(gridView, name, order, expr)
        elif gl == 2: # local function
            func = _localFunction(gridView, name, order, expr)
        elif gl == 3: # local function with self argument (i.e. from @gridFunction)
            func = _localFunction(gridView, name, order, lambda en,x: expr(en,x))
        return func
    # next we check for a cpp code snipett
    if isinstance(expr,str): # this is a cppFunction
        if gridView is None:
            raise AttributeError("for a function based on a C++ code a 'gridView' must be provided")
        if name=="_tmp" and fctName is None:
            raise AttributeError("the name of the C++ function has to be provided either through the 'name' or the 'fctName' parameter")
        if order is None:
            raise AttributeError("for a function based on a C++ code a 'order' must be provided")
        return cppFunction(gridView, name=name, order=order,
                     fctName=name,includes=io.StringIO(expr), **kwargs)

    # finally the decorator...
    if hasattr(expr,"hierarchicalGrid"):
        # this must be the gridView
        assert gridView is None, "don't provide also a gridView"
        gridView = expr
        expr = None
    elif expr is not None: # can't handle expr
        # can't handle this expression type
        raise AttributeError(f"can't handle given expr type {type(expr)}")
    if order is None:
        raise AttributeError("for a function based on Python code a 'order' must be provided")
    def gridFunction_decorator(func):
        from dune.fem.plotting import plotPointData
        gf = dune.grid.gridFunction(gridView,name=name,order=order)(func)
        setattr(gf.__class__,"as_ufl", lambda self: dune.ufl.GridFunction(self))
        setattr(gf.__class__,"integrate", lambda self: _uflFunction(gridView,self.name,order,self).integrate())
        # setattr(gf,"plot", lambda *args,**kwargs: plotPointData(gf,*args,**kwargs))
        gf.plot = lambda *args,**kwargs: plotPointData(gf,*args,**kwargs)
        return gf.as_ufl()
    return gridFunction_decorator

def levelFunction(gridView,name="levels"):
    @gridFunction(gridView,name=name,order=0)
    def levelFunction(e,x):
        return [e.level]
    return levelFunction

def partitionFunction(gridView,name="rank"):
    @dune.grid.GridFunction(gridView,name=name,order=0)
    class Partition(object):
        def __init__(self,rank):
            self.rank = rank
        def __call__(self,en,x):
            return [self.rank]
    return Partition(gridView.comm.rank)

def uflFunction(gridView, name, order, ufl, virtualize=True, scalar=False,
                predefined=None, *args, **kwargs):
    deprecated("dune.fem.function.uflFunction is deprecated. Use dune.fem.function.gridFunction(expr,gridView,name,order); note that the 'ufl' expression is first.")
    return _uflFunction(gridView,name,order,ufl,virtualize,scalar,predefined,*args,**kwargs)
def _uflFunction(gridView, name, order, ufl, virtualize=True, scalar=False,
                predefined=None, *args, **kwargs):
    expr = ufl
    import ufl
    if gridView is None:
        if type(expr) == list or type(expr) == tuple:
            expr = ufl.as_vector(expr)
        args, cc = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
        if len(args) > 0:
            raise AttributeError("no trial or test function should be included in the expression")
        gridView = set( c.gridView for c in cc if hasattr(c,"gridView") )
        if len(gridView) == 0:
            raise AttributeError("a 'gridView' has to be provided or the expression must contain a grid function.")
        if len(gridView) > 1:
            raise AttributeError("expression contains grid functions over different views.")
        gridView = gridView.pop() # there is only one element in this set - the unique gridView
    if order is None:
        try:
            order = ufl.algorithms.estimate_total_polynomial_degree(expr)
        except:
            raise AttributeError("a suitable order could not be determined - provide the order as argument")
    Func = dune.models.localfunction.UFLFunction(gridView, name, order,
            expr, renumbering=None,
            virtualize=virtualize,
            predefined=predefined, *args, **kwargs)
    if Func is None:
        raise AttributeError("could not generate ufl grid function from expression "+str(ufl))
    try:
        from dune.fem.plotting import plotPointData
        setattr(Func, "plot", plotPointData)
    except ImportError:
        setattr(Func, "plot", lambda *args,**kwargs:
           print("problem importing plotting utility - possibly matplotlib is missing?"))
    func = Func(gridView,name,order,*args,**kwargs)
    if not hasattr(func,"scalar"):
        func.scalar = scalar
    return func.as_ufl() if func is not None else None

def cppFunction(gridView, name, order, fctName, includes,
                *args,**kwargs):
    args = args + tuple(kwargs.get("args",()))
    gf = gridView.function(fctName,includes,*args,order=order,name=name)
    return dune.ufl.GridFunction( gf )

def _assertSizesMatch(space, dofVector):
    # only allow contiguous arrays as dof vectors
    try:
        assert dofVector.data.contiguous, "dofVector needs to provide a contiguous data memory, sliced arrays are not supported!"
    except AttributeError:
        pass

    # only allow arrays with one dimensional shapes
    try:
        assert len(dofVector.shape) == 1, "dofVector should be a simple array, i.e. len(shape) == 1"
    except AttributeError:
        pass

    # check that sizes and data size match (if vector size is larger always use the first part)
    try:
        assert space.size <= len(dofVector), f"space (size={space.size}) is larger than vector (size={len(dofVector)})!"
    except TypeError:
        assert space.size == dofVector.size, f"space (size={space.size}) is larger than vector (size={dofVector.size}) do not match!"

    if hasattr(dofVector, "dtype"):
        dtype = dofVector.dtype
        itemsize = np_dtype(dofVector.dtype).itemsize
        match = space._sizeOfField == itemsize and dtype == 'float64' if space.field == 'double' else True
        assert match, f"space (dtype={space._sizeOfField},{space.field}) and vector (dtype={itemsize},{dtype}) do not match!"

def discreteFunction(space, name, expr=None, dofVector=None):
    """create a discrete function

    Args:
        space: discrete function space
        name:  name of the discrete function
        expr:  analytical expression to interpolate

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    if dofVector is None:
        df = space.DiscreteFunction(space,name)
    else:
        _assertSizesMatch(space, dofVector)
        df = space.DiscreteFunction(name,space,dofVector)

    if expr is None and dofVector is None:
        df.clear()
    elif expr is not None:
        df.interpolate(expr)
    # df.scalar = space.scalar
    return df.as_ufl()

def tupleDiscreteFunction(*spaces, **kwargs):
    # from dune.fem.discretefunction import module, addAttr
    name = kwargs.get("name", "")
    try:
        tupleSpace = spaces[0]
        spaces = spaces[0].components
    except AttributeError:
        tupleSpace = dune.fem.space.tuple(*spaces)
    DiscreteFunction = tupleSpace.DiscreteFunction
    df = DiscreteFunction(tupleSpace,name)
    compNames = kwargs.get("components", None)
    if not compNames is None:
        components = df.components
        assert len(compNames) == len(components)
        for c, n in zip(components, compNames):
            setattr(df,n,c.as_ufl())
    return df.as_ufl()
