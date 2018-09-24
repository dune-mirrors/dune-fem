from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import warnings

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("GridAdaptation", "Dune::FemPy")

modules = {}

def module(grid):
    try:
        return modules[grid._typeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::GridAdaptation< " + grid._typeName + " >"
    includes = grid._includes + ["dune/fempy/py/grid/adaptation.hh"]
    moduleName = "adapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[grid._typeName] = module
    return module


def _adaptArguments(first,*args):
    try: # first see if first argument is a discrete function (should be only method)
        hgrid = first.grid.hierarchicalGrid
        args = list([*args,first])
    except AttributeError:
        pass
    try: # if its not a df test if it is a list/tuple
        if len(first)>1:
            assert len(args)==0,\
                   "only one list of discrete functions can be passed into the adaptation method"
            hgrid = first[0].grid.hierarchicalGrid
            args = first
    except TypeError: # okay apparently its a hgrid object (should be deprecated)
        hgrid = first
        if len(args)==1:
            args = args[0]
        warnings.warn("""
              passing in the hierarchical grid as first argument to the
              'dune.fem.adapt' function is not required and deprecated.
              """)
    # make sure all args are over the same grid
    assert all([a.grid.hierarchicalGrid==hgrid for a in args]),\
            "all discrete functions must be over the same hierarchical grid"
    return hgrid,args

def adapt(first, *args):
    hgrid,args = _adaptArguments(first,*args)
    module(hgrid).gridAdaptation(hgrid).adapt(args)

def loadBalance(first, *args):
    hgrid,args = _adaptArguments(first,*args)
    module(hgrid).gridAdaptation(hgrid).loadBalance(args)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
