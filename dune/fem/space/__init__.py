from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect
import sys

from dune.common.compatibility import isString
from dune.generator.generator import SimpleGenerator
from dune.fem import function

from ._spaces import *

def interpolate( self, func, name=None, **kwargs ):
    if not name:
        name = func.name
    return function.discreteFunction(self, name=name, expr=func, **kwargs)

def numpyfunction( self, data, name ):
    from dune.fem.function import numpyFunction
    return numpyFunction(self, data, name)

def storageToSolver(storage):
    if storage == "adaptive" or storage == "fem":
        return "fem"
    elif storage == "istl":
        return "istl"
    elif storage == "numpy":
        return "numpy"
    elif storage == "eigen":
        return "eigen"

generator = SimpleGenerator("Space", "Dune::FemPy")

def addAttr(module, cls, field, storage):
    setattr(cls, "_module", module)
    setattr(cls, "field", field )
    setattr(cls, "interpolate", interpolate )
    setattr(cls, "numpyfunction", numpyfunction )

    if not storage:
        storage = str("fem")
    if isString(storage):
        import dune.create as create
        storage = create.discretefunction( storageToSolver(storage) )(cls)
    else:
        storage = storage(cls)
    setattr(cls, "storage", storage )

fileBase = "femspace"

def module(field, storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Space, field, storage)
    return module
