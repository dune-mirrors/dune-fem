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

def interpolate(self, f):
    if ufl and isinstance(f, GridFunction):
        func = f.gf
    elif ufl and isinstance(f, ufl.core.expr.Expr):
        func = expression2GF(self.space.grid,f,self.space.order).as_ufl()
    else:
        try:
            gl = len(inspect.getargspec(f)[0])
            if gl == 1:   # global function
                func = function.globalFunction(self.space.grid, "tmp", self.space.order, f)
            elif gl == 2: # local function
                func = function.localFunction(self.space.grid, "tmp", self.space.order, f)
            elif gl == 3: # local function with self argument (i.e. from @gridFunction)
                func = function.localFunction(self.space.grid, "tmp", self.space.order, lambda en,x: f(en,x))
        except TypeError:
            func = f
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

fileBase = "femdiscretefunction"

def module(storage, includes, typeName, backend=None, *args, **kwargs):
    includes = includes + ["dune/fempy/py/discretefunction.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, bufferProtocol=True, **kwargs)
    addAttr(module, module.DiscreteFunction, storage)
    if not backend is None:
        addBackend(module.DiscreteFunction, backend)
    return module
