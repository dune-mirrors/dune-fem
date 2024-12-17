import hashlib
import inspect
import sys
import os
import functools
import numpy as np

from dune.deprecate import deprecated
import dune.common.module
from dune.common.utility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *
import ufl
from dune.ufl import GridFunction

def _uflToExpr(grid,order,f):
    if isinstance(f, list) or isinstance(f, tuple):
        if isinstance(f[0], ufl.core.expr.Expr):
            f = ufl.as_vector(f)
    if isinstance(f, GridFunction):
        return f
    if isinstance(f, ufl.core.expr.Expr):
        return function.gridFunction(f,gridView=grid,order=order).as_ufl()
    else:
        return f

def interpolate(space, expr, name=None, **kwargs):
    """interpolate a function into a discrete function space

    Args:
        space: discrete function space to interpolate into
        expr:  function to interpolate
        name:  name of the resulting discrete function

    Returns:
        DiscreteFunction: the constructed discrete function
    """
    expr = _uflToExpr(space.gridView,space.order,expr)
    if name is None:
        try:
            name = expr.name
        except AttributeError:
            raise ValueError("interpolation requires a name for the resulting discrete function - either the expression needs a name attribute or 'name' needs to be provided as second argument.")
    return function.discreteFunction(space, name=name, expr=expr, **kwargs)

def project(space, func, name=None, **kwargs):
    """project a (discrete) function into a discrete function space

    Args:
        space: discrete function space to project into
        func:  function to project
        name:  name of the resulting discrete function

    Returns:
        DiscreteFunction: the constructed discrete function
    """
    expr = _uflToExpr(space.gridView,space.order,func)
    if name is None:
        name = func.name
    # assert func.dimRange == space.dimRange, "range dimension mismatch"
    df = function.discreteFunction(space, name=name, expr=None, **kwargs)
    df.project(func)
    return df

def dfInterpolate(self, f):
    func = function.gridFunction(f, gridView=self.space.gridView, order=self.space.order).gf
    if func.dimRange != self.space.dimRange:
        raise AttributeError("trying to interpolate an expression"\
            " of size "+str(func.dimRange)+" into a space with range dimension = "\
            + str(self.space.dimRange))

    try:
        # API GridView: assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.dimRange == self.dimRange, "range dimension mismatch"
    except AttributeError:
        pass
    return self._interpolate(func)

def dfProject(self, f):
    func = function.gridFunction(f, gridView=self.space.gridView, order=self.space.order).gf
    if func.dimRange != self.space.dimRange:
        raise AttributeError("trying to interpolate an expression"\
            " of size "+str(func.dimRange)+" into a space with range dimension = "\
            + str(self.space.dimRange))

    try:
        # API GridView: assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.dimRange == self.dimRange, "range dimension mismatch"
    except AttributeError:
        pass
    return self._project(func)

def dfAssign(self, other):
    """ assign dofs of other discrete function to self
    Args:
        self: discrete function
        other: other discrete function
    """
    try:
        # try natural, when dof vectors match
        # (which is slightly more than df's match)
        self.dofVector.assign(other.dofVector)
    except:
        # this is the case when storage is different
        import io
        from dune.generator import algorithm
        _assignCode = \
"""
template <class A, class B>
void assignDF( A& a, const B& b )
{
  a.assign( b );
}
"""
        algorithm.run('assignDF', io.StringIO(_assignCode), self, other )

def localContribution(self, assembly):
    if assembly == "set":
        return self.setLocalContribution()
    elif assembly == "add":
        return self.addLocalContribution()
    else:
        raise ValueError("assembly can only be `set` or `add`")

def dummyPlot(*args,**kwargs):
    print("problem importing plotting utility - possibly matplotlib is missing?")

def backendProp(self):
    try:
        return self._backend
    except ImportError as ex:
        raise ex
    except:
        pass # still try numpy
    try:
        return np.array( self.dofVector, copy=False )
    except:
        pass
    return None

def addDFAttr(module, cls, backendName):
    cls._storage = property( lambda self: self._space.storage)
    cls.space = property( lambda self: self._space.as_ufl() )
    setattr(cls, "interpolate", dfInterpolate )
    setattr(cls, "assign", dfAssign )
    if hasattr(cls,"_project"):
        setattr(cls, "project", dfProject )
    setattr(cls, "localContribution", localContribution )
    cls.scalar = property(lambda self: self.space.scalar)
    try:
        from dune.fem.plotting import plotPointData
        setattr(cls, "plot", plotPointData)
    except ImportError:
        setattr(cls, "plot", dummyPlot)
    setattr(cls, backendName, property(backendProp))

# Issue: this property added later to the space will not be pickled
# so after loading the 'df' calling 'df.as_numpy' or similar will fail.
# Alternative is implemented above using the 'addDFAttr' during the
# DF class registration - those properties will be pickled.
# The problem is getting the passing in the correct backend name, e.g.,
# 'as_numpy' from the C++ code. A #define is used for that now but a better
# solution should be found in the future.
# Note: if the define is not set then the old approach is used.
def addBackend(Df,backend):
    if hasattr(Df,backend): # new approach has been used
        return
    def backend_(self):
        try:
            return self._backend
        except ImportError as ex:
            raise ex
        except:
            pass # still try numpy
        try:
            return np.array( self.dofVector, copy=False )
        except:
            pass
        return None
    setattr(Df,backend,property(backend_))

_defaultGenerator = SimpleGenerator(["Space","DiscreteFunction"], "Dune::FemPy")

def spcInterpolate(self,*args,**kwargs):
    return interpolate(self,*args,**kwargs)
def spcCodegen(self,moduleName,interiorQuadratureOrders, skeletonQuadratureOrders):
    return _codegen(self,moduleName,interiorQuadratureOrders, skeletonQuadratureOrders)
def spcProject(self,*args,**kwargs):
    return project(self,*args,**kwargs)
def spcFunction(self,*args,**kwargs):
    return function.discreteFunction(self,*args,**kwargs)
def zeroDF(space):
    try:
        if space._zero is not None:
            return space._zero
        raise AttributeError
    except AttributeError:
        space._zero = space.function(name="zero")
        space._zero.clear()
        return space._zero

def addAttr(module, self, field, scalar, codegen):
    setattr(self, "field", field)
    setattr(self, "scalar", scalar)
    setattr(self, "interpolate",functools.partial(spcInterpolate,self))
    setattr(self, "codegenStorage", codegen)
    setattr(self, "codegen", functools.partial(spcCodegen,self))
    module.Space.zero = property( zeroDF )

    DF = module.DiscreteFunction
    if hasattr(DF,"_project"):
        setattr(self, "project",functools.partial(spcProject,self))
    setattr(self, "function",functools.partial(spcFunction,self))

def addStorage(obj, storage):
    if not storage:
        storage = str("numpy")
    if isString(storage):
        import dune.create as create
        assert storage, "wrong storage (" + storage + ") passed to space"
        storage = create.discretefunction(storage)(obj)
    else:
        storage = storage(obj)
    setattr(obj, "storage", storage)
    return storage

fileBase = "femspace"

def addDiscreteFunction(space, storage):
    from dune.generator import Constructor

    storage = addStorage(space,storage)

    dfIncludes = storage.includes

    ctor = ()
    spaceType = space.cppTypeName
    if storage.name == "petsc":
        try:
            import petsc4py
            dfIncludes += [os.path.dirname(petsc4py.__file__)+"/include/petsc4py/petsc4py.h"]
            ctor = [Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle dofVector'],
                    ['if (import_petsc4py() != 0) {',
                     '  throw std::runtime_error("Error during import of petsc4py");',
                     '}',
                     'Vec petscVec = PyPetscVec_Get(dofVector.ptr());',
                     'typename DuneType::DofVectorType *dofStorage = new typename DuneType::DofVectorType(space,petscVec);',
                     'pybind11::cpp_function remove_petscStorage( [ dofStorage, dofVector, petscVec] ( pybind11::handle weakref ) {',
                     '  delete dofStorage;',
                     '  weakref.dec_ref();',
                     '} );',
                     'pybind11::weakref weakref( dofVector, remove_petscStorage );',
                     'weakref.release();',
                     'return new DuneType( name, space, *dofStorage );'],
                    ['"name"_a', '"space"_a', '"dofVector"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
                ]
        except:
            ctor = [Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle dofVector'],
                    ['std::cerr <<"Can not use constructor with dof vector argument because `petsc4py` was not found!\\n";',
                     'throw std::runtime_error("Can not use constructor with dof vector argument because `petsc4py` was not found!");',
                     'return new DuneType(name,space);'],
                    ['"name"_a', '"space"_a', '"dofVector"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
                   ]
    elif storage.name == "numpy" or storage.name == "fem":
        ctor = [Constructor(['const std::string &name', 'const ' +
            spaceType + '&space', 'pybind11::array_t<'+space.field+'> dofVector'],
                [space.field + ' *dof = static_cast< '+space.field+'* >( dofVector.request(false).ptr );',
                 'return new DuneType(name,space,dof);'],
                ['"name"_a', '"space"_a', '"dofVector"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
            ]
    elif storage.name == "istl":
        ctor = [Constructor(['const std::string &name', 'const ' +
            spaceType + '&space', 'typename DuneType::DofContainerType &dofVector'],
                ['return new DuneType(name,space,dofVector);'],
                ['"name"_a', '"space"_a', '"dofVector"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
            ]
    else:
        ctor = [Constructor(['const std::string &name', 'const ' + spaceType + '&space', 'pybind11::handle vec'],
                ['std::cerr <<"Can not use constructor with dof vector argument for this type of discrete function storage!\\n";',
                 'throw std::runtime_error("Can not use constructor with dof vector argument for this type of discrete function storage!");',
                 'return new DuneType(name,space);'],
                ['"name"_a', '"space"_a', '"vec"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
               ]
    dfIncludes += ["dune/fempy/py/discretefunction.hh"]
    return dfIncludes, storage.type, storage.backend, ctor

def addSpaceAttr(module,spc,field,scalar,codegen,storage,backend,clone):
    addAttr(module, spc, field, scalar, codegen)
    spc._backend = backend
    setattr(spc,"DiscreteFunction",module.DiscreteFunction)
    # add clone function if provided
    if clone is not None:
        setattr(spc,"clone",clone)
    addStorage(spc,storage)
    if not backend is None:
        addBackend(module.DiscreteFunction, backend)


def module(field, includes, typeName, *args,
           storage=None, scalar=False, codegen=False, clone=None,
           ctorArgs, generator=_defaultGenerator):

    includes = includes + ["dune/fempy/py/space.hh"]
    defines = []

    class DummySpace:
        cppTypeName = typeName
        cppIncludes = includes
    DummySpace.field = field
    dfIncludes, dfTypeName,  backend, dfArgs = addDiscreteFunction(DummySpace, storage)

    defines += [f"BACKENDNAME \"{backend}\""]

    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest() \
                          + "_" + hashlib.md5(dfTypeName.encode('utf-8')).hexdigest()

    module = generator.load(includes+dfIncludes, [typeName,dfTypeName], moduleName,
                            ((*args,), (*dfArgs,)),
                            dynamicAttr=[True,True],
                            bufferProtocol=[False,True],
                            options=[["std::shared_ptr<DuneType>"],[]],
                            defines=defines)

    spc = module.Space(*ctorArgs)
    addSpaceAttr(module,spc,field,scalar,codegen,storage,backend,clone)
    return spc

def _codegen(space,moduleName,interiorQuadratureOrders, skeletonQuadratureOrders):
    if interiorQuadratureOrders is None: interiorQuadratureOrders = []
    if skeletonQuadratureOrders is None: skeletonQuadratureOrders = []
    moduleNameExt = "i" + "".join(str(i) for i in interiorQuadratureOrders) + "_" +\
                    "s" + "".join(str(i) for i in skeletonQuadratureOrders)
    autogenFile = "autogeneratedcode_" + moduleName + "_" + moduleNameExt + ".hh"
    dune_py_dir   = dune.packagemetadata.getDunePyDir()
    generated_dir = dune_py_dir # os.path.join(dune_py_dir, 'python', 'dune', 'generated')
    codegenPath = generated_dir
    space._generateQuadratureCode(interiorQuadratureOrders,skeletonQuadratureOrders,codegenPath,autogenFile)
    return [autogenFile], moduleNameExt
