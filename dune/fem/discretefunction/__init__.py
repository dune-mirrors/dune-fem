from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import inspect
from dune.generator.generator import SimpleGenerator

from ._discretefunctions import *
from ._solvers import *

from dune.fem import function

try:
    import ufl
    from dune.ufl import GridFunction, expression2GF
except:
    pass

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
    else:
        try:
            if ufl and isinstance(func, GridFunction):
                func = func.gf
            elif ufl and isinstance(func, ufl.core.expr.Expr):
                func = expression2GF(self.space.grid,func,self.space.order)
        except NameError:
            pass
    return self._interpolate(func)

def localContribution(self, assembly):
    if assembly == "set":
        return self.setLocalContribution()
    elif assembly == "add":
        return self.addLocalContribution()
    else:
        raise ValueError("assembly can only be `set` or `add`")

def addAttr(module, cls, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)
    setattr(cls, "interpolate", interpolate )
    setattr(cls, "localContribution", localContribution )

fileBase = "femdiscretefunction"

def module(storage, includes, typeName, *args, **kwargs):
    includes = includes + ["dune/fempy/py/discretefunction.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, bufferProtocol=True, **kwargs)
    addAttr(module, module.DiscreteFunction, storage)
    return module
