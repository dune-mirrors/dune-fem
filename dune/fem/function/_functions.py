from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.models.localfunction

import dune.common.checkconfiguration as checkconfiguration

def registerGridFunctions(gridview):
    import hashlib
    from dune.generator import builder
    typeName = gridview._typeName
    moduleName = "femgridfunctions_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()

    includes = ["dune/fempy/py/grid/gridpart.hh", "dune/fempy/py/grid/function.hh"] + gridview._includes

    source = "".join(["#include <" + i + ">\n" for i in includes])
    source += "\n"
    source += "PYBIND11_PLUGIN( " + moduleName + " )\n"
    source += "{\n"
    source += "  typedef Dune::FemPy::GridPart< " + gridview._typeName + "> GridPart;\n"
    source += "  pybind11::module module( \"" + moduleName + "\" );\n"
    source += '  module.def( "globalGridFunction", Dune::FemPy::defGlobalGridFunction< GridPart >( module, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ));\n'
    source += '  module.def( "localGridFunction", Dune::FemPy::defLocalGridFunction< GridPart > ( module, "LocalGridFunction",  std::make_integer_sequence< int, 11 >() ));\n'
    source += "  return module.ptr();\n"
    source += "}\n"

    return builder.load(moduleName, source, "gridfunctions")

def globalFunction(gridview, name, order, value):
    module = registerGridFunctions(gridview)
    return module.globalGridFunction(gridview,name,order,value)

def localFunction(gridview, name, order, value):
    module = registerGridFunctions(gridview)
    return module.localGridFunction(gridview,name,order,value)

def levelFunction(gridview):
    return localFunction(gridview, "level", 0, lambda en,_: [en.level] )

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
    return dune.models.localfunction.UFLFunction(gridview, name, order, ufl, *args, **kwargs)

def discreteFunction(space, name, expr=None, *args, **kwargs):
    import dune.create as create
    df = create.discretefunction(space.storage,space,name,**kwargs)
    if expr:
        df.interpolate( expr )
    else:
        # df.interpolate( [0,]*df.dimRange )
        df.clear()
    return df

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
    includes = [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/fempy/py/common/numpyvector.hh" ] + space._module._includes
    spaceType = space._module._typeName
    field = space.field
    # typeName = "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
    typeName = "Dune::Fem::VectorDiscreteFunction< " +\
          spaceType + ", Dune::FemPy::NumPyVector< " + field + " > >"

    return module("numpy", includes, typeName).DiscreteFunction(space,name,vec)
