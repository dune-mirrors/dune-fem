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

def solve( scheme, target, rhs=None ):
    if rhs is None:
        if hasattr(scheme,"preconditioning") and scheme.preconditioning is not None:
            assert callable(scheme.preconditioning), "scheme.preconditioning needs to be a callable object: pre( u, v)!"
            return scheme._solve(target, scheme.preconditioning)
        else:
            return scheme._solve(target)
    else:
        return scheme._solve(rhs, target)


_defaultGenerator = SimpleGenerator("Scheme", "Dune::FemPy")

def addAttr(module, cls):
    setattr(cls, "solve", solve)

fileBase = "femscheme"

def module(includes, typeName, *args, backend=None,
           generator=_defaultGenerator):
    from dune.fem.space import addBackend
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True)
    addAttr(module, module.Scheme)
    JacobianOperator = getattr(module.Scheme,"JacobianOperator",None)
    if JacobianOperator is not None and hasattr(JacobianOperator,"_backend") and backend is not None:
        addBackend(JacobianOperator,backend)
    return module
