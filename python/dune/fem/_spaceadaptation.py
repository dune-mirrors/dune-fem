from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib

from dune.generator.generator import SimpleGenerator

_defaultGenerator = SimpleGenerator("SpaceAdaptation", "Dune::FemPy")

modules = {}

def module(space, generator=_defaultGenerator):
    storage, dfIncludes, dfTypeName, _, _, _ = space.storage
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
    dfs = {}
    for df in dfList:
        try:
            for dfc in df.components:
                try:
                    dfs[dfc.space] += [dfc]
                except KeyError:
                    dfs[dfc.space] = [dfc]
        except AttributeError:
            try:
                dfs[df.space] += [df]
            except KeyError:
                dfs[df.space] = [df]
    for s,df in dfs.items(): # if we really have different spaces then there is still the problem that the marker doesn't know which space we are on
        module(s).SpaceAdaptation(s).adapt(marker, df)
