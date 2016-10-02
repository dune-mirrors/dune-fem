from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect

from dune.generator.generator import SimpleGenerator
from dune.fem import discretefunction
from dune.fem.function import numpyFunction

from ._spaces import *

def interpolate( self, func, **kwargs ):
    import dune.create as create
    order = self.order
    try:
        gl = len(inspect.getargspec(func)[0])
    except:
        gl = 0
    if gl == 1:   # global function
        gf = create.function("global", self.grid, "tmp", order, func)
        return interpolate(self, gf, **kwargs)
    elif gl == 2: # local function
        gf = create.function("local", self.grid, "tmp", order, func)
        return interpolate(self, gf, **kwargs)
    elif gl == 0: # already a grid function
        storage = kwargs.pop('storage', "adaptive")
        try:
            name = kwargs.pop("name", func.name)
            df = create.discretefunction(storage, self, name=name, **kwargs)
        except:
            df = create.discretefunction(storage, self, **kwargs)
        df.interpolate(func)
        return df
    return None

def numpyfunction( self, data, name ):
    return numpyFunction(self, data, name)

generator = SimpleGenerator("Space", "Dune::FemPy")

def addAttr(module, cls, field):
    setattr(cls, "_module", module)
    setattr(cls, "field", field )
    setattr(cls, "interpolate", interpolate )
    setattr(cls, "numpyfunction", numpyfunction )

fileBase = "femspace"

def module(field, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/space.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Space, field)
    return module
