"""Functions for creating python modules and C++ classes for schemes.
"""

from __future__ import print_function
import importlib
import hashlib
import os.path
import re
# from termcolor import colored

from dune.deprecate import deprecated
from ._schemes import *

from dune.generator.generator import SimpleGenerator

@deprecated("the solve method on the 'scheme' now requires the 'target' argument, so first construct a discrete function using the 'interpolate' method on the 'space'")
def depSolve( scheme, target, rhs=None, name=None ):
    import dune.fem.function as function
    if target == None:
        if name == None:
            if hasattr(scheme, 'name'):
                name = scheme.name
            else:
                name = "solution"
        target = function.discreteFunction(scheme.space, name=name)
    if rhs is None:
        info = scheme._solve(target)
    else:
        info = scheme._solve(rhs, target)
    return target,info

def solve( scheme, target=None, rhs=None, name=None ):
    if target==None:
        return depSolve( scheme, target, rhs, name )
    else:
        if rhs is None:
            return scheme._solve(target)
        else:
            return scheme._solve(rhs, target)


generator = SimpleGenerator("Scheme", "Dune::FemPy")

def addAttr(module, cls):
    setattr(cls, "solve", solve)

fileBase = "femscheme"

def module(includes, typeName, *args, backend=None):
    from dune.fem.space import addBackend
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True)
    addAttr(module, module.Scheme)
    JacobianOperator = getattr(module.Scheme,"JacobianOperator",None)
    if JacobianOperator is not None and hasattr(JacobianOperator,"_backend") and backend is not None:
        addBackend(JacobianOperator,backend)
    return module
