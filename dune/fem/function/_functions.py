from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)


import dune.grid
import dune.fem.space
import dune.models.localfunction

import dune.common.checkconfiguration as checkconfiguration
from dune.common.hashit import hashIt

def registerGridFunctions(gridview):
    from dune.generator import builder
    typeName = gridview._typeName
    moduleName = "femgridfunctions_" + hashIt(typeName)

    includes = ["dune/fempy/py/grid/gridpart.hh", "dune/fempy/py/grid/function.hh"] + gridview._includes

    source = "#include <config.h>\n\n"
    source += "".join(["#include <" + i + ">\n" for i in includes])
    source += "\n"
    source += "PYBIND11_MODULE( " + moduleName + ", module )\n"
    source += "{\n"
    source += "  typedef Dune::FemPy::GridPart< " + gridview._typeName + "> GridPart;\n"
    source += '  module.def( "globalGridFunction", Dune::FemPy::defGlobalGridFunction< GridPart >( module, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ));\n'
    source += '  module.def( "localGridFunction", Dune::FemPy::defLocalGridFunction< GridPart > ( module, "LocalGridFunction",  std::make_integer_sequence< int, 11 >() ));\n'
    source += "}\n"

    return builder.load(moduleName, source, "gridfunctions")

def globalFunction(gridview, name, order, value):
    module = registerGridFunctions(gridview)
    return module.globalGridFunction(gridview,name,order,value).as_ufl()


def localFunction(gridview, name, order, value):
    module = registerGridFunctions(gridview)
    return module.localGridFunction(gridview,name,order,value).as_ufl()


def levelFunction(gridview):
    @dune.grid.gridFunction(gridview)
    def levelFunction(e,x):
        return [e.level]
    return levelFunction
    # return localFunction(gridview, "level", 0, lambda en,_: [en.level])


def partitionFunction(gridview):
    class Partition(object):
        def __init__(self,rank):
            self.rank = rank
        def __call__(self,en,x):
            return [self.rank]
    return localFunction(gridview, "rank", 0, Partition(gridview.comm.rank))


def cppFunction(gridview, name, order, code, *args, **kwargs):
    return dune.models.localfunction.generatedFunction(gridview, name, order, code, *args, **kwargs)


def uflFunction(gridview, name, order, ufl, *args, **kwargs):
    func = dune.models.localfunction.UFLFunction(gridview, name, order, ufl, *args, **kwargs)
    return func.as_ufl() if func is not None else None

def discreteFunction(space, name, expr=None, *args, **kwargs):
    """create a discrete function

    Args:
        space: discrete function space
        name:  name of the discrete function
        expr:  analytical expression to interpolate

    Returns:
        DiscreteFunction: the constructed discrete function
    """
    storage, dfIncludes, dfTypeName, _, _ = space.storage
    df = dune.fem.discretefunction.module(storage, dfIncludes, dfTypeName).DiscreteFunction(space,name)
    if expr is None:
        df.clear()
    else:
        df.interpolate(expr)
    return df.as_ufl()


def numpyFunction(space, vec, name="tmp", **unused):
    """create a discrete function - using the fem numpy storage as linear algebra backend
       Note: this is not a 'managed' discrete function, i.e., the storage
       is passed in and owned by the user. No resizing will take be done
       during grid modification.

    Args:
        space: discrete space
        vec: the vector storage (a numpy array)

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    from dune.fem.discretefunction import module
    assert vec.shape[0] == space.size, str(vec.shape[0]) +"!="+ str(space.size) + ": numpy vector has wrong shape"
    includes = [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/fempy/py/common/numpyvector.hh" ] + space._includes
    spaceType = space._typeName
    field = space.field
    typeName = "Dune::Fem::VectorDiscreteFunction< " +\
          spaceType + ", Dune::FemPy::NumPyVector< " + field + " > >"

    return module("numpy", includes, typeName).DiscreteFunction(space,name,vec).as_ufl()


def tupleDiscreteFunction(*spaces, **kwargs):
    from dune.fem.discretefunction import module, addAttr
    tupleSpace = dune.fem.space.combined(*spaces)
    dfIncludes = (space.storage[1] for space in spaces)
    dfTypeNames = (space.storage[2] for space in spaces)
    includes = sum(dfIncludes, ["dune/fem/function/tuplediscretefunction.hh"])
    typeName = "Dune::Fem::TupleDiscreteFunction< " + ", ".join(dfTypeNames) + " >"
    name = kwargs.get("name", "")
    df = module(tupleSpace.storage, includes, typeName, dynamicAttr=True).DiscreteFunction(tupleSpace, name)
    # create a discrete function for each space to ensure the DiscreteFunction is registered with pybind11
    for s in spaces:
        discreteFunction(s, "")
    #for i in range(df.components):
    #    addAttr(df._module, df[i].__class__, df._storage)
    compNames = kwargs.get("components", None)
    if not compNames is None:
        components = df.components
        assert len(compNames) == len(components)
        for c, n in zip(components, compNames):
            df.__dict__[n] = c
    return df.as_ufl()
