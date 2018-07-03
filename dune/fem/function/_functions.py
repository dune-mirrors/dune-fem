from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
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


def levelFunction(gridview,name="levels"):
    @dune.grid.gridFunction(gridview,name=name)
    def levelFunction(e,x):
        return [e.level]
    return levelFunction


def partitionFunction(gridview,name="rank"):
    class Partition(object):
        def __init__(self,rank):
            self.rank = rank
        def __call__(self,en,x):
            return [self.rank]
    return localFunction(gridview, name, 0, Partition(gridview.comm.rank))


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
    from dune.generator import Constructor
    storage, dfIncludes, dfTypeName, _, _,backend = space.storage

    if storage == "petsc":
        spaceType = space._typeName
        ctor = Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle vec'],
                ['if (import_petsc4py() != 0) {',
                 '  throw std::runtime_error("Error during import of petsc4py");',
                 '}',
                 'Vec petscVec = PyPetscVec_Get(vec.ptr());',
                 'typename DuneType::DofVectorType *dofStorage = new typename DuneType::DofVectorType(space,petscVec);',
                 '// std::cout << "setup_petscStorage " << dofStorage << " " << petscVec << std::endl;',
                 'pybind11::cpp_function remove_petscStorage( [ dofStorage, vec, petscVec] ( pybind11::handle weakref ) {',
                 '  // std::cout << "remove_petscStorage " << vec.ref_count() << " " << dofStorage << " " << petscVec << std::endl;',
                 '  delete dofStorage;',
                 '  weakref.dec_ref();',
                 '} );',
                 'pybind11::weakref weakref( vec, remove_petscStorage );',
                 'weakref.release();',
                 'return new DuneType( name, space, *dofStorage );'],
                ['"name"_a', '"space"_a', '"vec"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])
        DF = dune.fem.discretefunction.module(storage, dfIncludes, dfTypeName, backend, ctor)
        vec = kwargs.get("vec",None)
        if vec is None:
            df = DF.DiscreteFunction(space,name)
        else:
            df = DF.DiscreteFunction(name,space,vec)
            return df.as_ufl()
    else:
        df = dune.fem.discretefunction.module(storage, dfIncludes, dfTypeName, backend).DiscreteFunction(space,name)
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
    includes = space._includes + [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/fempy/py/common/numpyvector.hh" ]
    spaceType = space._typeName
    field = space.field
    typeName = "Dune::Fem::VectorDiscreteFunction< " +\
          spaceType + ", Dune::FemPy::NumPyVector< " + field + " > >"

    return module("numpy", includes, typeName, "as_numpy").DiscreteFunction(space,name,vec).as_ufl()

def petscFunction(space, vec, name="tmp", **unused):
    """create a discrete function - using the fem petsc storage as linear algebra backend
       Note: this is not a 'managed' discrete function, i.e., the storage
       is passed in and owned by the user. No resizing will take be done
       during grid modification.

    Args:
        space: discrete space
        vec: the vector storage (a numpy array)

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    return discreteFunction(space,name,vec=vec)
    from dune.generator import Constructor
    from dune.fem.discretefunction import module, petsc
    # assert vec.shape[0] == space.size, str(vec.shape[0]) +"!="+ str(space.size) + ": numpy vector has wrong shape"
    import petsc4py
    includes = space._includes +\
               [ os.path.dirname(petsc4py.__file__)+"/include/petsc4py/petsc4py.h" ] +\
               [ "dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh", "dune/fem/misc/petsc/petscvector.hh" ]

    spaceType = space._typeName
    typeName = "Dune::Fem::PetscDiscreteFunction< " + spaceType + ">"

    ctor = Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle vec'],
            ['std::cout << "UM: " << vec.ptr() << std::endl;',
             'if (import_petsc4py() != 0) {',
             '  std::cout << "ERROR: could not import petsc4py" << std::endl;',
             '  throw std::runtime_error("Error during import of petsc4py");',
             '}',
             'Vec petscVec = PyPetscVec_Get(vec.ptr());',
             'typename DuneType::DofVectorType *dofStorage = new typename DuneType::DofVectorType(space,petscVec);',
             'std::cout << "UM: " << petscVec << " " << *(dofStorage->getVector()) << std::endl;',
             'return new DuneType( name, space, *dofStorage );'],
            ['"name"_a', '"space"_a', '"vec"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])

    return module("petsc", includes, typeName, "as_petsc", ctor).DiscreteFunction(name,space,vec).as_ufl()


def tupleDiscreteFunction(*spaces, **kwargs):
    from dune.fem.discretefunction import module, addAttr
    try:
        tupleSpace = spaces[0]
        spaces = spaces[0].components
    except AttributeError:
        tupleSpace = dune.fem.space.tuple(*spaces)
    dfIncludes = (space.storage[1] for space in spaces)
    dfTypeNames = (space.storage[2] for space in spaces)
    includes = sum(dfIncludes, ["dune/fem/function/tuplediscretefunction.hh"])
    typeName = "Dune::Fem::TupleDiscreteFunction< " + ", ".join(dfTypeNames) + " >"
    name = kwargs.get("name", "")
    df = module(tupleSpace.storage, includes, typeName, dynamicAttr=True).DiscreteFunction(tupleSpace, name)
    # create a discrete function for each space to ensure the DiscreteFunction is registered with pybind11
    for s in spaces:
        discreteFunction(s, "")
    compNames = kwargs.get("components", None)
    if not compNames is None:
        components = df.components
        assert len(compNames) == len(components)
        for c, n in zip(components, compNames):
            df.__dict__[n] = c
    return df.as_ufl()
