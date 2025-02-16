from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib

from dune.generator.generator import SimpleGenerator
from dune.fem.deprecated import deprecated

_defaultGenerator = SimpleGenerator("SpaceAdaptation", "Dune::FemPy")

modules = {}

def module(space, generator=_defaultGenerator):

    dfTypeName = space.storage.type
    try:
        return modules[dfTypeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::SpaceAdaptation< " + space.storage.type + " >"
    includes = space.storage.includes + ["dune/fempy/py/spaceadapt.hh"]
    moduleName = "spaceadapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[dfTypeName] = module
    return module


def spaceAdapt(space, dfs, marker=None):
    """ Adapt the polynomial order according to the provided marker
        of the given discrete function space belonging to
        the discrete function or list of discrete functions. This function needs
        to be called for each discrete function space separately.

    Args:
        space: discrete function space to be adapted
        dfs:    discrete function or list of discrete functions belonging to space
        marker: marker function (same marker for all discrete functions).
                If marker is None it is assumed that marking was done already (by a user routine).

    Note: All discrete functions have to belong to the same space object.
          To adapt other discrete functions and spaces call spaceAdapt again.
          Otherwise different discrete function spaces have to be created.

    Returns:
        None
    """
    if not isinstance(dfs, (list,tuple)) and not hasattr(dfs, "space"):
        deprecated("spaceAdapt: old signature is used. Flip marker and discrete functions!")
        marker, dfs = dfs, marker

    if not isinstance(dfs, (list,tuple)):
        dfs = [dfs]

    # check that all dfs belong to the same space
    assert all([a.space==space for a in dfs]), "All discrete functions need to belong to the same discrete function space!"

    dfDict = {}
    for df in dfs:
        try:
            for dfc in df.components:
                try:
                    dfDict[dfc.space] += [dfc]
                except KeyError:
                    dfDict[dfc.space] = [dfc]
        except AttributeError:
            try:
                dfDict[df.space] += [df]
            except KeyError:
                dfDict[df.space] = [df]

    # adapt all dfs at once for one space
    # if we really have different spaces then there is still the problem that the marker doesn't know which space we are on
    if marker is None:
        for s,df in dfDict.items():
            module(s).SpaceAdaptation(s).adaptWithoutMarking(df)
    else:
        for s,df in dfDict.items():
            module(s).SpaceAdaptation(s).adapt(marker, df)
