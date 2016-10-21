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

def solve( scheme, target=None, name=None ):
    import dune.create as create
    if target == None:
        if name == None:
            if hasattr(scheme, 'name'):
                name = scheme.name
            else:
                name = "default"
        # target = create.discretefunction(scheme._storage, scheme.space, name=name)
        # target.interpolate( [0,]*scheme.dimRange )
        target = scheme.space.interpolate( lambda x:[0,]*scheme.dimRange, name=name)
    scheme._solve(target)
    return target

generator = SimpleGenerator("Scheme", "Dune::FemPy")

def storageToSolver(storage):
    if storage == "adaptive":
        return "fem"
    elif storage == "istl":
        return "istl"
    elif storage == "numpy":
        return "numpy"
    elif storage == "eigen":
        return "eigen"

def addAttr(module, cls):
    setattr(cls, "solve", solve)

fileBase = "femscheme"

def module(includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Scheme)
    return module
