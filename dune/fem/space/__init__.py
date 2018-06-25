from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect
import sys
import os

import dune.common.module
from dune.common.compatibility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *

import ufl
import dune.ufl


def interpolate(space, func, name=None, **kwargs):
    """interpolate a function into a discrete function space

    Args:
        space: discrete function space to interpolate into
        func:  function to interpolate
        name:  name of the resulting discrete function

    Returns:
        DiscreteFunction: the constructed discrete function
    """
    if name is None:
        name = func.name
    # assert func.dimRange == space.dimRange, "range dimension mismatch"
    return function.discreteFunction(space, name=name, expr=func, **kwargs)


def storageToSolver(storage):
    if storage == "adaptive":
        return "fem"
    else:
        return str(storage)

generator = SimpleGenerator("Space", "Dune::FemPy")

def addAttr(module, cls, field):
    setattr(cls, "_module", module)
    setattr(cls, "field", field)
    setattr(cls, "interpolate", interpolate)
    setattr(cls, "numpyFunction", function.numpyFunction)
    setattr(cls, "petscFunction", function.petscFunction)

def addStorage(obj, storage):
    if not storage:
        storage = str("fem")
    if isString(storage):
        import dune.create as create
        assert storageToSolver(storage), "wrong storage (" + storage + ") passed to space"
        storage = create.discretefunction(storageToSolver(storage))(obj)
    else:
        storage = storage(obj)
    setattr(obj, "storage", storage)

fileBase = "femspace"

def module(field, includes, typeName, *args,
           interiorQuadratureOrders=None, skeletonQuadratureOrders=None):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    defines = []
    if interiorQuadratureOrders is not None or\
       skeletonQuadratureOrders is not None:
        defines = ["USE_BASEFUNCTIONSET_CODEGEN"]
        includes = ["dune/fem/space/basisfunctionset/default_codegen.hh"] + includes
        moduleName = fileBase + "_" +\
            "i" + "".join(str(i) for i in interiorQuadratureOrders) + "_" +\
            "s" + "".join(str(i) for i in skeletonQuadratureOrders) + "_" +\
            hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True,
                            options=["std::shared_ptr<DuneType>"],
                            defines=defines)
    addAttr(module, module.Space, field)
    return module

def codegen(space,interiorQuadratureOrders, skeletonQuadratureOrders):
    if interiorQuadratureOrders is None: interiorQuadratureOrders = []
    if skeletonQuadratureOrders is None: skeletonQuadratureOrders = []
    dune_py_dir   = dune.common.module.get_dune_py_dir()
    generated_dir = dune_py_dir # os.path.join(dune_py_dir, 'python', 'dune', 'generated')
    codegenPath = generated_dir
    space._generateQuadratureCode(interiorQuadratureOrders,skeletonQuadratureOrders,codegenPath)
