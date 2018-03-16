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
    modules[dfTypeName] = module
    return module


def spaceAdapt(space, marker, dfList):
    if hasattr(space,"components"):
        # is a tuple space over a single space and all functions should be tuple dfs
        # singleSpace = space.components[0]
        dfs = {}
        for df in dfList:
            try:
                dfs[df.space] += [df]
            except KeyError:
                dfs[df.space] = [df]
        for s,df in dfs.items():
            print(s,df)
            module(s).SpaceAdaptation(s).adapt(marker, df)
    else:
        module(space).SpaceAdaptation(space).adapt(marker, dfList)
