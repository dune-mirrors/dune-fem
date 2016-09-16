"""Functions for creating python modules and C++ classes for schemes.
"""

from __future__ import print_function
import importlib
import subprocess
import hashlib
import os.path
import re
# from termcolor import colored

from dune.fem import comm
from .. import discretefunction

from dune.generator.generator import SimpleGenerator

def solve( scheme, rhs=None, target=None, name=None, assemble=True ):
    if name == None:
        if hasattr(scheme, 'name'):
            name = scheme.name
        else:
            name == "default"
    if target == None:
          target = discretefunction.create(scheme._storage, scheme.space, name=name)
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

def module(storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = "scheme_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.Scheme, storage)
    return module

schemeNames = { "h1"      : "dune.fem.scheme.h1",
                "dg"      : "dune.fem.scheme.dg",
                "mvdg"    : "dune.fem.scheme.nvdg",
                "stokes"  : "dune.fem.scheme.stokes",
                "burgers" : "dune.fen.scheme.burgers"
              }

def register(**kwargs):
    schemeNames.update(kwargs)

def create(scheme, *args, **kwargs):
    schemeModule = importlib.import_module(schemeNames[ scheme.lower() ])
    return schemeModule.create(*args,**kwargs)
