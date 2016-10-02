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

def solve( scheme, rhs=None, target=None, name=None, assemble=True ):
    import dune.create as create
    if name == None:
        if hasattr(scheme, 'name'):
            name = scheme.name
        else:
            name == "default"
    if target == None:
          target = create.discretefunction(scheme._storage, scheme.space, name=name)
          target.interpolate( [0,]*scheme.dimRange )
    if rhs == None:
        scheme._prepare()
    else:
        scheme._prepare(rhs)
    scheme._solve(target,assemble)
    return target

def spaceAndStorage(space_or_df,storage):
    try:
        space = space_or_df.space
    except:
        return space_or_df, storage
    assert storage == space_or_df._storage,\
               "missmatch in storages passed to scheme constructor -" +\
               "passed both a discrete function object and a storage" +\
               "argument that do not match: " +\
               storage + " != " + space_or_df._storage
               # colored( storage + " != " + space_or_df._storage, 'red')
    storage = space_or_df._storage
    return space, storage

generator = SimpleGenerator("Scheme", "Dune::FemPy")

def addAttr(module, cls, storage):
    setattr(cls, "solve", solve)
    setattr(cls, "_storage", storage)

fileBase = "femscheme"

def module(storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Scheme, storage)
    return module
