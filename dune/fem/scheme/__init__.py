"""Functions for creating python modules and C++ classes for schemes.
"""

from __future__ import print_function
import importlib
import subprocess
import hashlib
import os.path
import re
# from termcolor import colored

from ._schemes import *

from dune.generator.generator import SimpleGenerator

def solve( scheme, rhs=None, target=None, name=None ):
    import dune.fem.function as function
    if target == None:
        if name == None:
            if hasattr(scheme, 'name'):
                name = scheme.name
            else:
                name = "default"
        target = function.discreteFunction(scheme.space, name=name)
    if rhs is None:
        info = scheme._solve(target)
    else:
        info = scheme._solve(rhs, target)
    return target,info

generator = SimpleGenerator("Scheme", "Dune::FemPy")

def addAttr(module, cls):
    setattr(cls, "solve", solve)

fileBase = "femscheme"

def module(includes, typeName, *args):
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True)
    addAttr(module, module.Scheme)
    return module
