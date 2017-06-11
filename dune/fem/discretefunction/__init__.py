from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import inspect
from dune.generator.generator import SimpleGenerator

from ._discretefunctions import *
from ._solvers import *

from dune.fem import function

generator = SimpleGenerator("DiscreteFunction", "Dune::FemPy")

def interpolate(self, func):
    try:
        gl = len(inspect.getargspec(func)[0])
    except TypeError:
        gl = 0
    if gl == 1:   # global function
        func = function.globalFunction(self.space.grid, "tmp", self.space.order, func)
    elif gl == 2: # local function
        func = function.localFunction(self.space.grid, "tmp", self.space.order, func)
    return self._interpolate(func.gf)

def addAttr(module, cls, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)
    setattr(cls, "interpolate", interpolate )

fileBase = "femdiscretefunction"

def module(storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/discretefunction.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods, bufferProtocol=True)
    addAttr(module, module.DiscreteFunction, storage)
    return module
