from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("DiscreteFunction", "Dune::FemPy")

def addAttr(module, cls, storage):
    setattr(cls, "_module", module)
    setattr(cls, "_storage", storage)

def module(storage, includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/discretefunction.hh"]
    moduleName = "discretefunction_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.DiscreteFunction, storage)
    return module

discretefunctionNames = { "adaptive" : "dune.fem.discretefunction.adaptive",
                          "fem"      : "dune.fem.discretefunction.adaptive",
                          "istl"     : "dune.fem.discretefunction.istl",
                          "numpy"    : "dune.fem.discretefunction.numpy",
                        }

def register(**kwargs):
    discretefunctionNames.update(kwargs)

def create(discretefunction, *args, **kwargs):
    discretefunctionModule = importlib.import_module(discretefunctionNames[ discretefunction.lower() ])
    return discretefunctionModule.create(*args,**kwargs)

#############################################

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
