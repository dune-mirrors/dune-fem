import hashlib
import inspect
import sys
import os

from dune.deprecate import deprecated
import dune.common.module
from dune.common.utility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *

try:
    import ufl
    from dune.ufl import GridFunction, expression2GF
except:
    pass

def _uflToExpr(grid,order,f):
    if not ufl: return f
    if isinstance(f, list) or isinstance(f, tuple):
        if isinstance(f[0], ufl.core.expr.Expr):
            f = ufl.as_vector(f)
    if isinstance(f, GridFunction):
        return f
    if isinstance(f, ufl.core.expr.Expr):
        return expression2GF(grid,f,order).as_ufl()
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
    if isinstance(f, list) or isinstance(f, tuple):
        if isinstance(f[0], ufl.core.expr.Expr):
            f = ufl.as_vector(f)
    dimExpr = 0
    if isinstance(f, GridFunction):
        func = f.gf
        dimExpr = func.dimRange
    elif isinstance(f, ufl.core.expr.Expr):
        func = expression2GF(self.space.gridView,f,self.space.order).as_ufl()
        if func.ufl_shape == ():
            dimExpr = 1
        else:
            dimExpr = func.ufl_shape[0]
    else:
        try:
            gl = len(inspect.getargspec(f)[0])
            func = None
        except TypeError:
            func = f
            if isinstance(func,int) or isinstance(func,float):
                dimExpr = 1
            else:
                dimExpr = len(func)
        if func is None:
            if gl == 1:   # global function
                func = function.globalFunction(self.space.gridView, "tmp", self.space.order, f).gf
            elif gl == 2: # local function
                func = function.localFunction(self.space.gridView, "tmp", self.space.order, f).gf
            elif gl == 3: # local function with self argument (i.e. from @gridFunction)
                func = function.localFunction(self.space.gridView, "tmp", self.space.order, lambda en,x: f(en,x)).gf
            dimExpr = func.dimRange

    if dimExpr == 0:
        raise AttributeError("can not determine if expression shape"\
                " fits the space's range dimension")
    elif dimExpr != self.space.dimRange:
        raise AttributeError("trying to interpolate an expression"\
            " of size "+str(dimExpr)+" into a space with range dimension = "\
            + str(self.space.dimRange))

    try:
        # API GridView: assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.gridView == self.gridView, "can only interpolate with same grid views"
        assert func.dimRange == self.dimRange, "range dimension mismatch"
    except AttributeError:
        pass
    return self._interpolate(func)

def dfProject(self, f):
    if isinstance(f, list) or isinstance(f, tuple):
        if isinstance(f[0], ufl.core.expr.Expr):
            f = ufl.as_vector(f)
    dimExpr = 0
    if isinstance(f, GridFunction):
        func = f.gf
        dimExpr = func.dimRange
    elif isinstance(f, ufl.core.expr.Expr):
        func = expression2GF(self.space.gridView,f,self.space.order).as_ufl()
        if func.ufl_shape == ():
            dimExpr = 1
        else:
            dimExpr = func.ufl_shape[0]
    else:
        try:
            gl = len(inspect.getargspec(f)[0])
            func = None
        except TypeError:
            func = f
            if isinstance(func,int) or isinstance(func,float):
                dimExpr = 1
            else:
                dimExpr = len(func)
        if func is None:
            if gl == 1:   # global function
                func = function.globalFunction(self.space.gridView, "tmp", self.space.order, f).gf
            elif gl == 2: # local function
                func = function.localFunction(self.space.gridView, "tmp", self.space.order, f).gf
            elif gl == 3: # local function with self argument (i.e. from @gridFunction)
                func = function.localFunction(self.space.gridView, "tmp", self.space.order, lambda en,x: f(en,x)).gf
            dimExpr = func.dimRange

    if dimExpr == 0:
        raise AttributeError("can not determine if expression shape"\
                " fits the space's range dimension")
    elif dimExpr != self.space.dimRange:
        raise AttributeError("trying to interpolate an expression"\
            " of size "+str(dimExpr)+" into a space with range dimension = "\
            + str(self.space.dimRange))
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

def addDFAttr(module, cls, spc, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)
    cls.space = property( lambda self: self._space.as_ufl() )
    setattr(cls, "interpolate", dfInterpolate )
    setattr(cls, "assign", dfAssign )
    if hasattr(cls,"_project"):
        setattr(cls, "project", dfProject )
    setattr(cls, "localContribution", localContribution )
    try:
        from dune.fem.plotting import plotPointData
        setattr(cls, "plot", plotPointData)
    except ImportError:
        setattr(cls, "plot", lambda *args,**kwargs:
           print("problem importing plotting utility - possibly matplotlib is missing?"))

def addBackend(Df,backend):
    def backend_(self):
        try:
            return self._backend
        except:
            pass
        try:
            import numpy as np
            return np.array( self.dofVector, copy=False )
        except:
            pass
        return None
    setattr(Df,backend,property(backend_))

_defaultGenerator = SimpleGenerator(["Space","DiscreteFunction"], "Dune::FemPy")

def addAttr(module, self, field, scalar, codegen):
    setattr(self, "_module", module)
    setattr(self, "field", field)
    setattr(self, "scalar", scalar)
    setattr(self, "codegenStorage", codegen)
    setattr(self, "interpolate",
            lambda *args,**kwargs: interpolate(self,*args,**kwargs))
    setattr(self, "codegen",
            lambda moduleName,interiorQuadratureOrders, skeletonQuadratureOrders:
                  _codegen(self,moduleName,interiorQuadratureOrders, skeletonQuadratureOrders)
           )

    DF = module.DiscreteFunction
    if hasattr(DF,"_project"):
        setattr(self, "project",
            lambda *args,**kwargs: project(self,*args,**kwargs))
    setattr(self, "function",
            lambda *args,**kwargs: function.discreteFunction(self,*args,**kwargs))
    DF.scalar = property(lambda self: self.space.scalar)

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
    storage, dfIncludes, dfTypeName, _, _,backend = addStorage(space,storage)

    ctor = ()
    spaceType = space.cppTypeName
    if storage == "petsc":
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
    elif storage == "numpy" or storage == "fem":
        ctor = [Constructor(['const std::string &name', 'const ' +
            spaceType + '&space', 'pybind11::array_t<double> dofVector'],
                ['double *dof = static_cast< double* >( dofVector.request(false).ptr );',
                 'return new DuneType(name,space,dof);'],
                ['"name"_a', '"space"_a', '"dofVector"_a', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()'])
            ]
    elif storage == "istl":
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
    return dfIncludes, dfTypeName, backend, ctor

def module(field, includes, typeName, *args,
           storage=None, scalar=False, codegen=False, ctorArgs,
           generator=_defaultGenerator):
    includes = includes + ["dune/fempy/py/space.hh"]
    defines = []

    class DummySpace:
        cppTypeName = typeName
        cppIncludes = includes
    DummySpace.field = field
    dfIncludes, dfTypeName,  backend, dfArgs = addDiscreteFunction(DummySpace, storage)

    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest() \
                          + "_" + hashlib.md5(dfTypeName.encode('utf-8')).hexdigest()

    module = generator.load(includes+dfIncludes, [typeName,dfTypeName], moduleName,
                            ((*args,), (*dfArgs,)),
                            dynamicAttr=[True,True],
                            bufferProtocol=[False,True],
                            options=[["std::shared_ptr<DuneType>"],[]],
                            defines=defines)

    spc = module.Space(*ctorArgs)
    addAttr(module, spc, field, scalar, codegen)
    setattr(spc,"DiscreteFunction",module.DiscreteFunction)
    addDFAttr(module, module.DiscreteFunction, spc, addStorage(spc,storage))
    if not backend is None:
        addBackend(module.DiscreteFunction, backend)
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
