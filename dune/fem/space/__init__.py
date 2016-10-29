from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect

from dune.generator.generator import SimpleGenerator
from dune.fem import discretefunction,function

from ._spaces import *

def interpolate( self, func, name=None, **kwargs ):
    import dune.create as create
    order = self.order
    try:
        gl = len(inspect.getargspec(func)[0])
    except:
        gl = 0
    if gl == 1:   # global function
        gf = function.globalFunction(self.grid, "tmp", order, func)
        return interpolate(self, gf, name, **kwargs)
    elif gl == 2: # local function
        gf = function.localFunction(self.grid, "tmp", order, func)
        return interpolate(self, gf, name, **kwargs)
    elif gl == 0: # already a grid function
        storage = self.storage
        if not name:
            name = func.name
        df = function.discreteFunction(self, name=name, expr=func, **kwargs)
        return df
    return None

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
    if isinstance(storage,str):
        import dune.create as create
        storage = create.discretefunction( storageToSolver(storage) )(cls)
    setattr(cls, "storage", storage )

fileBase = "femspace"

def module(field, storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Space, field, storage)
    return module
