from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect
import sys

from dune.common.compatibility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *


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
    return function.discreteFunction(space, name=name, expr=func, **kwargs)


def storageToSolver(storage):
    if storage == "adaptive" or storage == "fem":
        return "fem"
    elif storage == "istl":
        return "istl"
    elif storage == "numpy":
        return "numpy"
    elif storage == "eigen":
        return "eigen"
    elif storage == "petsc":
        return "petsc"

generator = SimpleGenerator("Space", "Dune::FemPy")

def addAttr(module, cls, field, storage):
    setattr(cls, "_module", module)
    setattr(cls, "field", field)
    setattr(cls, "interpolate", interpolate)
    setattr(cls, "numpyFunction", function.numpyFunction)

    if not storage:
        storage = str("fem")
    if isString(storage):
        import dune.create as create
        assert storageToSolver(storage), "wrong storage (" + storage + ") passed to space"
        storage = create.discretefunction(storageToSolver(storage))(cls)
    else:
        storage = storage(cls)
    setattr(cls, "storage", storage)

fileBase = "femspace"

def module(field, storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Space, field, storage)
    return module
