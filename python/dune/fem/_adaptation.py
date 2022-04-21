from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib, warnings

try:
    import ufl
    from dune.ufl import GridFunction, expression2GF
    from ufl.algorithms.analysis import extract_arguments_and_coefficients
except ImportError:
    ufl = None

from dune.generator.generator import SimpleGenerator

_defaultGenerator = SimpleGenerator("GridAdaptation", "Dune::FemPy")

modules = {}

import logging, traceback
logger = logging.getLogger(__name__)

def module(grid, generator=_defaultGenerator):
    try:
        return modules[grid.cppTypeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::GridAdaptation< " + grid.cppTypeName + " >"
    includes = grid.cppIncludes + ["dune/fempy/py/grid/adaptation.hh"]
    moduleName = "adapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[grid.cppTypeName] = module
    return module


def _adaptArguments(first,*args):
    try: # first see if first argument is a discrete function (should be only method)
        hgrid = first.space.gridView.hierarchicalGrid
        args = list([*args,first])
        return hgrid,args
    except AttributeError:
        pass
    if isinstance(first,list) or isinstance(first,tuple):
        assert len(args)==0,\
               "only one list of discrete functions can be passed into the adaptation method"
        hgrid = first[0].space.gridView.hierarchicalGrid
        args = first
    else: # okay apparently its a hgrid object (should be deprecated if args is not empty)
        hgrid = first
        assert hasattr(hgrid,"leafView"), "second argument is neither a discrete function nor a hierarchical grid"
        if len(args)==0:
            return hgrid,()
        if len(args)==1:
            args = args[0]
        warnings.warn("""
              passing in the hierarchical grid as first argument to the
              'dune.fem.adapt' function is not required and deprecated.
              """)
    return hgrid,args

def adapt(first, *args):
    hgrid,args = _adaptArguments(first,*args)

    # make sure all args are over the same grid
    assert all([a.space.gridView.hierarchicalGrid==hgrid for a in args]),\
            "all discrete functions must be over the same hierarchical grid"
    # make sure all gridview can be adapted
    try:
        adapt  = all([a.space.gridView.canAdapt==True for a in args])
        # adapt &= not any([a.space.storage[0]=="petsc" for a in args])
    except AttributeError:
        adapt = False
    assert adapt,\
            "the grid views for all discrete functions need to support adaptivity e.g. `adpative` view also `petsc' storage can not be used at this point in time"

    module(hgrid).gridAdaptation(hgrid).adapt(args)

def loadBalance(first, *args):
    hgrid,args = _adaptArguments(first,*args)
    module(hgrid).gridAdaptation(hgrid).loadBalance(args)

def mark(indicator, refineTolerance, coarsenTolerance=0.,
         minLevel=0, maxLevel=None, minVolume=-1.0, maxVolume=-1.0,
         gridView=None,
         markNeighbors=False):
    if ufl and (isinstance(indicator, list) or isinstance(indicator, tuple)):
        indicator = ufl.as_vector(indicator)
    if ufl and isinstance(indicator, ufl.core.expr.Expr):#
        if gridView is None:
            _, coeff_ = extract_arguments_and_coefficients(indicator)
            gridView = [c.grid for c in coeff_ if hasattr(c,"grid")]
            gridView += [c.gridView for c in coeff_ if hasattr(c,"gridView")]
            if len(gridView) == 0:
                raise ValueError("if a ufl expression is passed as indicator then the 'gridView' must also be provided")
            gridView = gridView[0]
        indicator = expression2GF(gridView,indicator,0)
    if gridView is None:
        gridView = indicator.gridView
    try:
        if not gridView.canAdapt:
            raise AttributeError("indicator function must be over grid view that supports adaptation")
    except AttributeError:
        raise AttributeError("indicator function must be over grid view that supports adaptation")

    if maxLevel==None:
        maxLevel = -1
    return gridView.mark(indicator,refineTolerance,coarsenTolerance,minLevel,maxLevel, minVolume, maxVolume, markNeighbors)

def markNeighbors(indicator, refineTolerance, coarsenTolerance=0.0,
                  minLevel=0, maxLevel=None,
                  minVolume=-1.0, maxVolume=-1.0,
                  gridView=None):
    return mark(indicator, refineTolerance,
                coarsenTolerance=coarsenTolerance,
                minLevel=minLevel,
                maxLevel=maxLevel,
                minVolume=minVolume,
                maxVolume=maxVolume,
                gridView=gridView,
                markNeighbors=True)

def doerflerMark(indicator, theta, maxLevel=None, layered=0.05):
    try:
        if not indicator.space.gridView.canAdapt:
            raise AttributeError("indicator function must be over grid view that supports adaptation")
    except:
        raise AttributeError("indicator function must be over grid view that supports adaptation")
    if maxLevel==None:
        maxLevel = -1
    return indicator.space.gridView.doerflerMark(indicator,theta,maxLevel,layered)

def globalRefine(level, first, *args):
    hgrid,args = _adaptArguments(first,*args)

    # make sure all args are over the same grid
    assert all([a.gridView.hierarchicalGrid==hgrid for a in args]),\
            "all discrete functions must be over the same hierarchical grid"
    # make sure all gridview can be adapted
    try:
        adapt  = all([a.gridView.canAdapt==True for a in args])
    except AttributeError:
        adapt = False
    assert adapt,\
            "the grid views for all discrete functions need to support adaptivity e.g. `adpative` view"

    module(hgrid).gridAdaptation(hgrid).globalRefine(level, args)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
