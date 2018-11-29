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
    raise NotImplementedError("cpp function is not working at the moment")
    return dune.models.localfunction.generatedFunction(gridview, name, order, code, *args, **kwargs)


def uflFunction(gridview, name, order, ufl, virtualize=True, *args, **kwargs):
    Func = dune.models.localfunction.UFLFunction(gridview, name, order,
            ufl, renumbering=None,
            virtualize=virtualize, *args, **kwargs)
    if Func is None:
        raise AttributeError("could not generate ufl grid function from expression "+str(ufl))
    func = Func(gridview,name,order,*args,**kwargs)
    return func.as_ufl() if func is not None else None

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
        df = space.DiscreteFunction(name,space,dofVector)

    if expr is None and dofVector is None:
        df.clear()
    elif expr is not None:
        df.interpolate(expr)
    return df.as_ufl()

    from dune.generator import Constructor
    storage, dfIncludes, dfTypeName, _, _,backend = space.storage

    if storage == "petsc":
        spaceType = space._typeName
        try:
            import petsc4py
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
                    ['"name"_a', '"space"_a', '"vec"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
        except:
            ctor = Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle vec'],
                    ['std::cerr <<"Can not use constructor with dof vector argument because `petsc4py` was not found!\\n";',
                     'throw std::runtime_error("Can not use constructor with dof vector argument because `petsc4py` was not found!");',
                     'return new DuneType(name,space);'],
                    ['"name"_a', '"space"_a', '"vec"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
        DF = dune.fem.discretefunction.module(storage, dfIncludes, dfTypeName, backend, ctor)
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
