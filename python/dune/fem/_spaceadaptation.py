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

class SpaceMarker:
    def __init__(self, indicator, refineTolerance, coarsenTolerance=0.,
                 minOrder = 0,
                 maxOrder = -1,
                 markNeighbors = False,
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
        self._setIndicator(indicator)
        self.refineTolerance  = refineTolerance
        self.coarsenTolerance = coarsenTolerance
        self.minOrder         = minOrder
        self.maxOrder         = maxOrder
        self.markNeighbors    = markNeighbors
        self.statistics       = statistics

    def _getIndicator(self):
        return self._indicator

    def _setIndicator(self, indicator):
        self._indicator = indicator
        # if indicator changes reset gridFunction
        self._indicatorGF = None

    # indicator attribute
    indicator = property(fget=_getIndicator, fset=_setIndicator)

    def __call__(self, space):
        """ Perform marking of given space

        Args:
            space: discrete function space to be marked

        Returns:
            tuple( int, int ) if statistics was set to True.

        """

        assert space.canAdapt, "Provided space does not support p-adaptation!"
        # turn indicator into gridFunction if not already gridFunction
        if not self._indicatorGF:
            self._indicatorGF = _indicatorToGridFunction( self.indicator, space.gridView )

        assert self.minOrder >= 0 and self.minOrder <= space.order, "Invalid value for minOrder, must be between 0 and space.order, including!"
        if self.maxOrder == -1: self.maxOrder = space.order
        assert self.maxOrder >= self.minOrder and self.maxOrder <= space.order, "Invalid value for maxOrder, must be between minOrder and space.order, including!"

        return space._mark(self._indicatorGF, self.refineTolerance, self.coarsenTolerance,
                           self.minOrder, self.maxOrder,
                           self.markNeighbors, self.statistics )


def spaceAdapt(marker, dfs, *args, **kwargs):
    """ Adapt the polynomial order of the discrete function space according to the provided marker (or a
        previously set marking). The discrete functions provided will be adjusted accordingly.

    Args:
        marker: marker function (same marker for all discrete functions).
        dfs:    discrete function or list of discrete functions belonging to space

    Note: To adapt discrete functions with different markers call this routine
          multiple times with different lists of discrete functions.

    Returns:
        None.
    """

    from dune.ufl import FemSpace
    if isinstance(marker, FemSpace):
        deprecated("spaceAdapt: old signature used. Do not pass space anymore, use spaceAdapt( marker, dfs )!")
        assert len(args) == 1, "Signature failure"
        return spaceAdapt( dfs, args[0] )

    if not isinstance(dfs, (list,tuple)):
        dfs = [dfs]

    # untangle all different spaces for the list of dfs
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
    if isinstance(marker, SpaceMarker):
        def _adapt(spc, dfList):
            # perform marking
            marker( spc )
            # adapt space and dfs
            module(spc).SpaceAdaptation(spc).adapt(dfList)
    else:
        def _adapt( spc, dfList ):
            assert s.canAdapt, "Provided space or sub-space is not p-adaptive!"
            module(s).SpaceAdaptation(s).adapt(marker, dfList)

    # perform adaptation
    for s,dfList in dfDict.items():
        _adapt( s, dfList );
