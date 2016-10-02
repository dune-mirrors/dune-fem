from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import inspect

from dune.generator.generator import SimpleGenerator
from dune.fem import discretefunction
from dune.fem.discretefunction import numpy as numpydf

from ._spaces import *

def interpolate( self, func, **kwargs ):
    from dune.create import create
    order = self.order
    try:
        gl = len(inspect.getargspec(func)[0])
    except:
        gl = 0
    if gl == 1:   # global function
        return interpolate(self,self.grid.globalGridFunction("gf", order, func), **kwargs)
    elif gl == 2: # local function
        return interpolate(self,self.grid.localGridFunction("gf", order, func), **kwargs)
    elif gl == 0: # already a grid function
        storage = kwargs.pop('storage', "Adaptive")
        try:
            name = kwargs.pop("name", func.name)
            df = create.discretefunction(storage, self, name=name, **kwargs)
        except:
            df = create.discretefunction(storage, self, **kwargs)
        df.interpolate(func)
        return df
    return None

def numpyfunction( self, name, data ):
    assert data.shape[0] == self.size(), "vector has wrong shape"
    return numpydf(self, data, name)

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
