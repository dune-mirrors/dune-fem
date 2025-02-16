from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib

from dune.generator.generator import SimpleGenerator
from dune.fem.deprecated import deprecated
from dune.fem._adaptation import _indicatorToGridFunction

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

def _spaceMark(space, indicator, refineTolerance, coarsenTolerance=0.,
               minOrder=0, maxOrder=None,
               markNeighbors=False,
               statistics = False ):
    """
    Mark the polynomial order on each element for either increase or decrease or
    keep. The function space.adapt has to be called subsequently.

    Args:
        indicator         A piecewise constant grid function (or ufl expression) holding the estimated error (or indicator value)
        refineTolerance   An element E will be marked for order increase if indicator(E) > refineTolerance.
        coarsenTolerance  An element E will be marked for order decrease if indicator(E) < coarsenTolerance.
        minOrder          Minimal order to be assumed (in [0, space.order).
        maxOrder          Maximal order to be assumed (in [minOrder, space.order]).
        markNeighbors     If True, for an element E that is marked for order increase all it's neighbors will be marked for at least the same order.
        statistics        If True, the returned tuple contains the number of elements marked for increase and decrease, respectively.

    Returns:
        tuple( int, int )

    """
    indicator = _indicatorToGridFunction( indicator, space.gridView )

    assert minOrder >= 0 and minOrder <= space.order, "Invalid value for minOrder, must be between 0 and space.order, including!"

    if maxOrder is None: maxOrder = space.order
    assert maxOrder >= minOrder and maxOrder <= space.order, "Invalid value for maxOrder, must be between minOrder and space.order, including!"

    return space._mark(indicator, refineTolerance, coarsenTolerance, minOrder,
                       maxOrder, markNeighbors, statistics)


def _spaceAdapt(space, dfs, marker=None):
    """ Adapt the polynomial order of the discrete function space according to the provided marker (or a
        previously set marking). The discrete functions provided will be adjusted accordingly.

    Args:
        dfs:    discrete function or list of discrete functions belonging to space
        marker: marker function (same marker for all discrete functions).
                If marker is None it is assumed that marking was done already (by a user routine).

    Note: All discrete functions have to belong to this space object.
          To adapt other discrete functions and spaces call adapt on that
          object.

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
            assert s.canAdapt, "Provided space or sub-space is not p-adaptive!"
            module(s).SpaceAdaptation(s).adaptWithoutMarking(df)
    else:
        for s,df in dfDict.items():
            assert s.canAdapt, "Provided space or sub-space is not p-adaptive!"
            module(s).SpaceAdaptation(s).adapt(marker, df)

def spaceAdapt(space, *args, **kwargs):
    deprecated("spaceAdapt: deprecated method, use space.adapt directly!")
    space.adapt( *args, **kwargs )
