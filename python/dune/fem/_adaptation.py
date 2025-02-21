from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib, warnings

import ufl
from dune.fem.function import gridFunction
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from dune.generator.generator import SimpleGenerator
from dune.fem.deprecated import deprecated

_defaultGenerator = SimpleGenerator("GridAdaptation", "Dune::FemPy")

modules = {}

import logging
logger = logging.getLogger(__name__)

def module(grid, generator=_defaultGenerator):
    try:
        return modules[grid.cppTypeName]
    except KeyError:
        pass

    typeName = "Dune::FemPy::GridAdaptation< " + grid.cppTypeName + " >"
    includes = grid.cppIncludes + ["dune/fempy/py/grid/adaptation.hh"]
    moduleName = "adapt_" + hashlib.md5(typeName.encode('utf8')).hexdigest()

    module = generator.load(includes, typeName, moduleName)
    modules[grid.cppTypeName] = module
    return module


class AdaptationMarkerBase:
    @staticmethod
    def indicatorToGridFunction( indicator, gridView=None ):
        if isinstance(indicator, (list, tuple)):
            indicator = ufl.as_vector(indicator)
        if isinstance(indicator, ufl.core.expr.Expr):
            if gridView is None:
                _, coeff_ = extract_arguments_and_coefficients(indicator)
                gridView = [c.grid for c in coeff_ if hasattr(c,"grid")]
                gridView += [c.gridView for c in coeff_ if hasattr(c,"gridView")]
                if len(gridView) == 0:
                    raise ValueError("if a ufl expression is passed as indicator then the 'gridView' must also be provided")
                gridView = gridView[0]
            indicator = gridFunction(indicator,gridView=gridView,order=0)
        return indicator

    # constructor
    def __init__(self, indicator, refineTolerance, coarsenTolerance=0.,
                 markNeighbors = False, statistics = False):
        self._setIndicator(indicator)
        self.refineTolerance  = refineTolerance
        self.coarsenTolerance = coarsenTolerance
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

    def indicatorGF(self, gridView = None):
        if not self._indicatorGF:
            self._indicatorGF = AdaptationMarkerBase.indicatorToGridFunction( self.indicator, gridView )
        return self._indicatorGF

    def __call__(self, adaptee):
        raise Exception("needs to be overloaded by derived class!")
### end AdaptationMarkerBase

class GridMarker(AdaptationMarkerBase):
    def __init__(self, indicator, refineTolerance, coarsenTolerance=0.,
                 minLevel=0, maxLevel=-1,
                 minVolume=-1.0, maxVolume=-1.0,
                 gridView=None,
                 markNeighbors=False,
                 statistics = False,
                 strategy='default',
                 layered = 0.05):
        """
        Mark elements of the underlying hierarchical grid for refinement or
        coarsening. The function fem.adapt has to be called subsequently.

        Args:
            indicator         A piecewise constant grid function (or ufl expression) holding the estimated error (or indicator value)
            refineTolerance   An element E will be marked for refinement if indicator(E) > refineTolerance.
            coarsenTolerance  An element E will be marked for coarsening if indicator(E) < coarsenTolerance.
            minLevel          An element E will only be marked for coarsening if E.level > minLevel.
            maxLevel          An element E will only be marked for refinement if E.level < maxLevel.
            minVolume         An element E will only be marked for refinement if E.geometry.volume > minVolume.
            maxVolume         An element E will only be marked for coarsening if E.geometry.volume < maxVolume.
            gridView          gridView (and hierarchical grid) if it cannot be extracted from indicator.
            markNeighbors     If True, for an element E that is marked for refinement all it's neighbors will be marked for refinement also.
            statistics        If True, the __call__ operator returns a tuple that contains the number of elements marked for refinement and coarsening, respectively.
            strategy          Selected strategy, available are 'default' and 'doerfler'.
            layered           Parameter for Doerfler layered strategy. 0 selects standard Doerfler strategy.
        """

        super().__init__(indicator, refineTolerance, coarsenTolerance,
                         markNeighbors, statistics)

        self.minLevel  = minLevel
        if maxLevel == None: maxLevel = -1
        self.maxLevel  = maxLevel
        self.minVolume = minVolume
        self.maxVolume = maxVolume
        self.gridView  = gridView
        self.strategy  = strategy # default, doerfler
        self.layered   = layered

    def __call__(self):
        # get indicator grid function
        indicatorGF = self.indicatorGF(self.gridView)
        # obtain grid view
        if self.gridView is None:
            self.gridView = indicatorGF.gridView
        try:
            if not self.gridView.canAdapt:
                raise AttributeError("indicator function must be over grid view that supports adaptation")
        except AttributeError:
            raise AttributeError("indicator function must be over grid view that supports adaptation")

        return self.gridView.mark(indicatorGF, self.refineTolerance, self.coarsenTolerance,
                                  self.minLevel, self.maxLevel, self.minVolume, self.maxVolume, self.markNeighbors, self.statistics)
### end GridMarker

def _adaptArguments(first,*args):
    try: # first see if first argument is a discrete function (should be only method)
        hgrid = first.space.gridView.hierarchicalGrid
        args = list([*args,first])
        return hgrid,args
    except AttributeError:
        pass
    if isinstance(first,list) or isinstance(first,tuple):
        assert len(args)==0,\
               "only one list of discrete functions can be passed into the adaptation method"
        hgrid = first[0].space.gridView.hierarchicalGrid
        args = first
    else: # okay apparently its a hgrid object (should be deprecated if args is not empty)
        hgrid = first
        assert hasattr(hgrid,"leafView"), "second argument is neither a discrete function nor a hierarchical grid"
        if len(args)==0:
            return hgrid,()
        if len(args)==1:
            args = args[0]
        warnings.warn("""
              passing in the hierarchical grid as first argument to the
              'dune.fem.adapt' function is not required and deprecated.
              """)

    # make sure all args are over the same grid
    assert all([a.space.gridView.hierarchicalGrid==hgrid for a in args]),\
            "all discrete functions must be over the same hierarchical grid"
    # make sure all gridview can be adapted
    try:
        adapt  = all([a.space.gridView.canAdapt==True for a in args])
        # adapt &= not any([a.space.storage[0]=="petsc" for a in args])
    except AttributeError:
        adapt = False
    assert adapt, "the grid views for all discrete functions need to support adaptivity"
    return hgrid,args

def gridAdapt(marker, *args):
    """ Adapt the underlying hierarchical grid of the discrete function passed
        as arguments.

    Args:
        marker: A marking callable to mark the underlying hierarchical grid.
        *args: a single discrete function or a list or tuple of
        discrete functions which should be projected to the new grid.

    Note: All discrete functions have to belong to the same hierarchical grid.

    Returns:
        None
    """
    from dune.ufl import GridFunction
    # check if marker is discrete function and if so call with marker=None
    if marker is not None and (not callable(marker) or isinstance(marker, GridFunction)):
        return adapt(None, marker, *args)

    assert len(args) >= 1
    first = args[0]
    args  = args[1:]

    # obtain hierarchical grid and discrete functions
    hgrid,args = _adaptArguments(first,*args)

    # mark grid entities if marker is not None, otherwise just adapt grid
    if marker is not None:
        if isinstance(marker, GridMarker):
            # GridMarker object using gridView.mark
            marker()
        else:
            # Python or C++ callable
            hgrid.mark( marker )

    # perform adaptation step
    module(hgrid).gridAdaptation(hgrid).adapt(args)

def adapt(marker, *args):
    # deprecated("adapt: call gridAdapt instead!")
    return gridAdapt( marker, *args )


def loadBalance(first, *args):
    """ Balance the underlying hierarchical grid such that the work load is as
        good as possible distributed to the available MPI processes.

    Args:
        (first, *args): a single discrete function or a list or tuple of
        discrete functions which should be projected to the new grid.

    Note: All discrete functions have to belong to the same hierarchical grid.

    Returns:
        None
    """
    hgrid, args = _adaptArguments(first,*args)
    module(hgrid).gridAdaptation(hgrid).loadBalance(args)

def mark(indicator, refineTolerance, coarsenTolerance=0.,
         minLevel=0, maxLevel=None,
         minVolume=-1.0, maxVolume=-1.0,
         gridView=None,
         markNeighbors=False,
         statistics = False ):
    """
    Mark elements of the underlying hierarchical grid for refinement or
    coarsening. The function fem.adapt has to be called subsequently.

    Args:
        indicator         A piecewise constant grid function (or ufl expression) holding the estimated error (or indicator value)
        refineTolerance   An element E will be marked for refinement if indicator(E) > refineTolerance.
        coarsenTolerance  An element E will be marked for coarsening if indicator(E) < coarsenTolerance.
        minLevel          An element E will only be marked for coarsening if E.level > minLevel.
        maxLevel          An element E will only be marked for refinement if E.level < maxLevel.
        minVolume         An element E will only be marked for refinement if E.geometry.volume > minVolume.
        maxVolume         An element E will only be marked for coarsening if E.geometry.volume < maxVolume.
        gridView          gridView (and hierarchical grid) if it cannot be extracted from indicator.
        markNeighbors     If True, for an element E that is marked for refinement all it's neighbors will be marked for refinement also.
        statistics        If True, the returned tuple contains the number of elements marked for refinement and coarsening, respectively.

    Returns:
        tuple( int, int )

    """
    marker = GridMarker( indicator, refineTolerance, coarsenTolerance=coarsenTolerance,
                         minLevel=minLevel, maxLevel=maxLevel, minVolume=minVolume, maxVolume=maxVolume,
                         gridView=gridView, markNeighbors=markNeighbors, statistics=statistics)
    return marker()

def markNeighbors(indicator, refineTolerance, coarsenTolerance=0.0,
                  minLevel=0, maxLevel=None,
                  minVolume=-1.0, maxVolume=-1.0,
                  gridView=None):
    return mark(indicator, refineTolerance,
                coarsenTolerance=coarsenTolerance,
                minLevel=minLevel,
                maxLevel=maxLevel,
                minVolume=minVolume,
                maxVolume=maxVolume,
                gridView=gridView,
                markNeighbors=True)

def doerflerMark(indicator, theta, maxLevel=None, layered=0.05):
    """
    Mark elements of the underlying hierarchical grid for refinement or
    coarsening based on the Doerfler strategy. The function fem.adapt has to be called subsequently.
    See also: W. Dörfler, A Convergent Adaptive Algorithm for Poisson's Equation,
    SIAM J. Numer. Anal. 33 (3), 1106-1124, 1996.

    Args:
        indicator         A piecewise constant grid function holding the estimated error (or indicator value)
        theta             Tolerance for sum estimator value.
        maxLevel          An element E will only be marked for refinement if E.level < maxLevel (default is None).
        layered           Parameter for layered Doerfler marking (default is 0.05).

    Returns:
        tuple( int, int )

    """
    try:
        if not indicator.space.gridView.canAdapt:
            raise AttributeError("indicator function must be over grid view that supports adaptation")
    except:
        raise AttributeError("indicator function must be over grid view that supports adaptation")
    if maxLevel==None:
        maxLevel = -1
    return indicator.space.gridView.doerflerMark(indicator,theta,maxLevel,layered)

def globalRefine(level, first, *args):
    """ Adapt the underlying hierarchical grid uniformly by the level and adjust
        the discrete function passed as arguments.

    Args:
        level: Number of global refinements
        (first, *args): a single discrete function or a list or tuple of
        discrete functions which should be projected to the new grid.

    Note: All discrete functions have to belong to the same hierarchical grid.

    Returns:
        None
    """
    hgrid,args = _adaptArguments(first,*args)

    # if no discrete functions were supplied then simply globally refine grid
    if len(args) == 0:
        hgrid.globalRefine(level)
        return
    module(hgrid).gridAdaptation(hgrid).globalRefine(level, args)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
