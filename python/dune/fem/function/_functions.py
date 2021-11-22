from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
import logging
logger = logging.getLogger(__name__)
from numpy import dtype as np_dtype

import dune.ufl
import dune.grid
import dune.fem.space
import dune.models.localfunction
from dune.generator import builder

import dune.common.checkconfiguration as checkconfiguration
from dune.common.hashit import hashIt

def integrate(grid,expression,order):
    try:
        return expression.integrate()
    except AttributeError:
        return uflFunction(grid,"tmp",order,expression).integrate()
    # return dune.ufl.expression2GF(grid,expression,order).integrate()
# perhaps a general
#    def assemble(grid/space, expression, order):
# would be better. If expression is a function this would return a
# fieldvector with the integral of this function, if expression is a
# functional a discreteFunctional and for a bilinear form a matrix is
# returned. This could also include a "piecewise" option that (at least in
# the case of a function) returns the integral on each element.


def globalFunction(gridView, name, order, value):
    gf = gridView.function(value,name=name,order=order)
    return dune.ufl.GridFunction(gf)
def localFunction(gridView, name, order, value):
    gf = gridView.function(value,name=name,order=order)
    return dune.ufl.GridFunction(gf)
# decorators similar to dune.python but with ufl support
def gridFunction(view,name,order):
    def gridFunction_decorator(func):
        gf = dune.grid.gridFunction(view,name=name,order=order)(func)
        setattr(gf.__class__,"as_ufl", lambda self: dune.ufl.GridFunction(self))
        setattr(gf.__class__,"integrate", lambda self: uflFunction(view,"tmp",order,self).integrate())
        return gf.as_ufl()
    return gridFunction_decorator
# this is not going to work - needs fixing
# def GridFunction(view, name=None):
#     def GridFunction_decorator(cls):
#         return dune.ufl.GridFunction(dune.grid.GridFunction(view,name,order)(cls))
#     return GridFunction_decorator


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
    Func = dune.models.localfunction.UFLFunction(gridView, name, order,
            ufl, renumbering=None,
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

    # check that sizes and data size match
    try:
        assert space.size == len(dofVector), f"space (size={space.size}) and vector (size={len(dofVector)}) do not match!"
    except TypeError:
        assert space.size == dofVector.size, f"space (size={space.size}) and vector (size={dofVector.size}) do not match!"

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
            setattr(df,n,c)
    return df.as_ufl()
