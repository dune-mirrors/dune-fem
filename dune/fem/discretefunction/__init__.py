from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
from dune.generator.generator import SimpleGenerator

from ._discretefunctions import *
from ._solvers import *

generator = SimpleGenerator("DiscreteFunction", "Dune::FemPy")

def addAttr(module, cls, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)

fileBase = "femdiscretefunction"

def module(storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/discretefunction.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.DiscreteFunction, storage)
    return module
