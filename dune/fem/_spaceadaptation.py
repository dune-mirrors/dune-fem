from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("SpaceAdaptation", "Dune::FemPy")

modules = {}

def module(space):
    storage, dfIncludes, dfTypeName, _, _ = space.storage
    try:
        return modules[dfTypeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::SpaceAdaptation< " + dfTypeName + " >"
    includes = dfIncludes + ["dune/fempy/py/spaceadapt.hh"]
    moduleName = "spaceadapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[grid._typeName] = module
    return module


def spaceAdapt(space, marker, *args, **kwargs):
    module(space).spaceAdaptation(space).adapt(*args, **kwargs)
