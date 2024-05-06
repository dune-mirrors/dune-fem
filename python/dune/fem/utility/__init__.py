from __future__ import print_function
from dune.generator import algorithm, path
import io
from dune.common import comm

class EmptyResult:
    def __init__(self, obj):
        self.obj = obj
    def result(self):
        return self.obj

# empty class making code work that uses ThreadPoolExecutor
class EmmptyThreadPoolExecutor:
    def __init__(self,*args, **kwargs):
        pass

    # return instance of this class upon entering the with statement
    def __enter__(self, *args, **kwargs):
        return EmmptyThreadPoolExecutor()

    # do nothing here
    def __exit__(self, *args, **kwargs):
        pass

    # return an object holding the result
    def submit(self, function, *args, **kwargs):
        return EmptyResult( function(*args, **kwargs) )

# in serial runs use ThreadPoolExecutor, otherwise the empty dummy
if comm.size == 1:
    from concurrent.futures import ThreadPoolExecutor as FemThreadPoolExecutor
else:
    # for parallel runs this seems not to work, so avoid it
    FemThreadPoolExecutor = EmmptyThreadPoolExecutor

###################################################################
###################################################################
##
## GridWidth: Convenience class to compute grid width 'h'.
##
###################################################################
###################################################################
class GridWidth:
    def __init__(self, gridView):
        self._gridView = gridView
        code = """
#include <memory>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/gridpart/common/gridpartadapter.hh>

template <class GridView>
std::pair<double, int > gridWidth(const GridView& gv, const double h, const int sequence)
{
  typedef typename GridView::Grid GridType;
  typedef Dune::Fem::DofManager<GridType> DofManagerType;

  Dune::Fem::GridPartAdapter< GridView > gp( gv );
  const int currentSequence = DofManagerType :: instance( gp.grid() ).sequence();
  if( currentSequence != sequence )
  {
    double newWidth = Dune::Fem::GridWidth::calcGridWidth(gp);
    return std::make_pair( newWidth, currentSequence );
  }
  else
  {
    return std::make_pair( h, sequence );
  }
}
"""
        self._h = 0.0
        self._sequence = -1 # default to trigger re-computation
        self._gridWidth = algorithm.load("gridWidth", io.StringIO(code), self._gridView, self._h, self._sequence )

    def __call__(self):
        """
        Returns:
        --------
            h computed as min( h_E ) forall E in gridView where h_E = |E| / min(|e|) for all e in E
        """
        # potentially update h if grid has changed
        self._h, self._sequence = self._gridWidth(self._gridView, self._h, self._sequence )
        return self._h

def gridWidth( gridView ):
    """
    Computes the minimal grid width which over all elements by computing min { |E|/|e| }
    for all intersections e of an element.

    Parameters:
        gridView  grid view to compute gridwidth h for. h is comp

    Returns: h computed as min(|E|/|e|) forall E in gridView
    """
    return GridWidth(gridView)()

###################################################################
###################################################################
##
## inspectBoundaryIds: Project boundary ids of a given grid view to
##                     a piecewise constant function
##
###################################################################
###################################################################
def inspectBoundaryIds( gridView ):
    code = """
#include <dune/fem/misc/boundaryidprovider.hh>

template <class DF>
void inspectBoundary( DF& df )
{
  Dune::Fem::projectBoundaryIds( df );
}
"""
    from dune.fem.space import finiteVolume
    space = finiteVolume( gridView )
    bndIds = space.interpolate([0], name="bndId" )
    algorithm.run("inspectBoundary", io.StringIO(code), bndIds )
    return bndIds

###################################################################
###################################################################
##
## Sampler: Line, Point, and Boundary sample.
##
###################################################################
###################################################################
class Sampler:
    def __init__(self, gridFunction ):
        self.gridFunction = gridFunction
        self.lineSampler = None
        self.pointSampler = None
        self.boundarySampler = None

    def lineSample(self,x0,x1,N, allowOutside=False):
        from dune.common import FieldVector
        import numpy
        x0, x1 = FieldVector(x0), FieldVector(x1)
        if self.lineSampler is None:
            _lineSampleCode = \
            """
            #ifndef FEMPY_UTILITY_HH
            #define FEMPY_UTILITY_HH
            #include <vector>
            #include <utility>
            #include <dune/fem/misc/linesegmentsampler.hh>
            #include <dune/fem/gridpart/common/entitysearch.hh>
            #include <dune/fem/function/localfunction/const.hh>

            template <class GF, class DT>
            std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
            sample(const GF &gf, DT &start, DT &end, int n, bool allowOutside)
            {
              Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
              std::vector<DT> coords(n);
              std::vector<typename GF::RangeType> values(n);
              if (allowOutside)
                sampler(gf,values, std::nothrow);
              else
                sampler(gf,values);
              sampler.samplePoints(coords);
              return std::make_pair(coords,values);
            }
            #endif

            """
            self.lineSampler = algorithm.load('sample', io.StringIO(_lineSampleCode),
                                 self.gridFunction, x0, x1, N, allowOutside)

        p,v = self.lineSampler( self.gridFunction, x0, x1, N, allowOutside )
        if self.gridFunction.scalar:
            x,y = numpy.zeros(len(p)), numpy.zeros(len(p))
        else:
            x,y = numpy.zeros(len(p)), numpy.zeros(( len(p),self.gridFunction.dimRange) )
        length = (x1-x0).two_norm
        for i in range(len(x)):
            x[i] = (p[i]-x0).two_norm / length
            if self.gridFunction.scalar:
                y[i] = v[i][0]
            else:
                y[i] = v[i]
        return x,y

    def pointSample(self, x0):
        from dune.common import FieldVector
        import numpy
        x0 = FieldVector(x0)
        if self.pointSampler is None:
            if hasattr(self.gridFunction,"gridView"):
                _pointSampleCode = \
                """
                #ifndef FEMPY_UTILITY_HH
                #define FEMPY_UTILITY_HH
                #include <dune/fem/misc/linesegmentsampler.hh>
                #include <dune/fem/gridpart/common/entitysearch.hh>
                #include <dune/fem/function/localfunction/const.hh>

                template <class GF, class DT>
                typename GF::RangeType sample(const GF &gf, DT &point)
                {
                  typedef typename GF::DiscreteFunctionSpaceType::GridPartType GridPartType;
                  Dune::Fem::EntitySearch<GridPartType> search(gf.space().gridPart());
                  const auto &entity = search(point);
                  const auto localPoint = entity.geometry().local(point);
                  return constLocalFunction(gf,entity).evaluate(localPoint);
                }
                #endif // FEMPY_UTILITY_HH

                """
                self.pointSampler = algorithm.load('sample', io.StringIO(_pointSampleCode), self.gridFunction, x0)
            else:
                _pointSampleCode = \
                """
                #ifndef FEMPY_UTILITY_HH
                #define FEMPY_UTILITY_HH
                #include <tuple>
                #include <dune/fem/gridpart/common/entitysearch.hh>
                #include <dune/fempy/py/grid/gridpart.hh>

                template <class GV, class DT>
                auto sample(const GV &gv, DT &point)
                {
                  const auto& gp = Dune::FemPy::gridPart(gv);
                  Dune::Fem::EntitySearch<Dune::FemPy::GridPart<GV>>
                    search(gp);
                  const auto &entity = search(point);
                  const auto localPoint = entity.geometry().local(point);
                  return std::make_tuple(entity,localPoint);
                }
                #endif // FEMPY_UTILITY_HH

                """
                self.pointSampler = algorithm.load('sample', io.StringIO(_pointSampleCode), self.gridFunction, x0)

        v = self.pointSampler(self.gridFunction, x0 )
        try:
            if self.gridFunction.scalar:
                return v[0]
            else:
                return v
        except AttributeError:
            return v

    def boundarySample(self, boundaryId=None, boundaryDomain=None, order=-1):
        """
        Parameter:
           boundaryId      integer or list or tuple of integer identifier of the boundary segments to be sampled,
                           can be none if a boundaryDomain was provided.
           boundaryDomain  characteristic function `boundaryDomain( x, boudnaryId ) -> Boolean` of the boundary domain to be sampled,
                           can be none if a boundary id was provided.
           order           polynomial order of the quadrature to be used (default
                           is the grid functions polynomial order.

        Return:
           points, values  numpy arays holding the coordinates and evaluation of the grid function.
                           Note, if the boundary segment provided it not found,
                           points and values will be empty.

        """
        import numpy

        if boundaryId is None and boundaryDomain is None:
            raise ValueError("Either a boundaryId or a boundaryDomain has to be provided in order to identify the boundary domain to be sampled.")
        # set default boundary id to zero, then it will be ignored
        if boundaryId is None:
            boundaryId = [-1]

        # integers will be converted to list with one element
        if isinstance(boundaryId,int):
            boundaryId = [boundaryId]

        # make sure boundary ids are either a list or tuple structure
        # otherwise the passing to C++ will not work
        assert isinstance(boundaryId, (list,tuple))

        # set default boundary domain when not given (selects all points passed)
        hasBoundaryDomain = boundaryDomain is not None
        if not hasBoundaryDomain:
            def bnd(x, bndId):
                return True
            boundaryDomain = bnd

        boundaryDomain.cppIncludes = ["dune/common/fvector.hh"]
        boundaryDomain.cppTypeName = "std::function<bool(const Dune::FieldVector<double, " + str(self.gridFunction.gridView.dimensionworld) + ">&)>"

        if self.boundarySampler is None:
            _bndSampleCode = \
"""
#include <algorithm>
#include <vector>
#include <utility>
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>

template <class GF>
std::pair<std::vector<typename GF::DomainType>, std::vector<typename GF::RangeType>>
sampleBnd(const GF &gf, const std::vector<int>& boundaryIds,
          const std::function<bool(const typename GF::DomainType&)>& boundaryDomain,
          const bool hasBoundaryDomain, const int optionalOrder)
{
  typedef typename GF::DomainType DomainType;
  typedef typename GF::RangeType  RangeType;

  std::vector<DomainType> coords;
  std::vector<RangeType> values;

  std::vector<int> ret(1, 0);

  const int order = (optionalOrder == -1) ? gf.order() : optionalOrder;

  auto pointIsValid = [&ret, &boundaryDomain, &hasBoundaryDomain](const typename GF::DomainType& x)
  {
    /*
    if( hasBoundaryDomain )
    {
      boundaryDomain( x, ret );
      return ret[0] == 1;
    }
    else
    */
      return true;
  };

  typedef typename GF::GridPartType GridPartType;
  const GridPartType& gridPart = gf.gridPart();

  const auto boundaryIdsEnd = boundaryIds.end();
  Dune::Fem::ConstLocalFunction< GF > localGF( gf );
  for( const auto &element : Dune::elements( gridPart, Dune::Partitions::interior ) )
  {
    if( element.hasBoundaryIntersections() )
    {
      for( const auto &intersection : Dune::intersections( gridPart, element ) )
      {
        if( !intersection.boundary() ) continue;

        const int bndId = Dune::Fem::boundaryId< GridPartType >(intersection);
        const bool hasBoundaryId = (std::find(boundaryIds.begin(), boundaryIdsEnd, bndId ) !=  boundaryIdsEnd);
        if( hasBoundaryId || hasBoundaryDomain )
        {
          auto guard = Dune::Fem::bindGuard( localGF, element );
          auto geometry = element.geometry();
          // sample points
          typedef  Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
          FaceQuadratureType quad(gridPart,intersection, order, FaceQuadratureType::INSIDE);
          const int nop = quad.nop();

          coords.reserve( coords.size() + nop );
          values.reserve( values.size() + nop );

          RangeType value;
          for( int qp=0; qp<nop; ++qp )
          {
            const auto point = geometry.global(Dune::Fem::coordinate( quad[ qp ] ) );
            if( pointIsValid( point ) )
            {
              coords.push_back( point );
              localGF.evaluate( quad[ qp ], value );
              values.push_back( value );
            }
          }
        }
      }
    }
  }

  return std::make_pair(coords,values);
}
"""
            self.boundarySampler = algorithm.load('sampleBnd', io.StringIO(_bndSampleCode), self.gridFunction, boundaryId, boundaryDomain, hasBoundaryDomain, order )

        p,v = self.boundarySampler( self.gridFunction, boundaryId, boundaryDomain, hasBoundaryDomain, order )
        # points could be empty if boundaryId was not found

        if len(p) > 0 and len(p[0]) > 1:
            x = numpy.zeros( (len(p), len(p[0]) ) )
        else:
            x = numpy.zeros( len(p) )

        if self.gridFunction.scalar:
            y = numpy.zeros(len(p))
        else:
            y = numpy.zeros((len(p),self.gridFunction.dimRange) )

        for i in range(len(x)):
            x[i] = p[i]
            if self.gridFunction.scalar:
                y[i] = v[i][0]
            else:
                y[i] = v[i]
        return x,y

#########################################################################
# end Sampler
#########################################################################


def lineSample(gridFunction,x0,x1,N, allowOutside=False):
    return Sampler(gridFunction).lineSample(x0, x1, N, allowOutside)

def pointSample(gridFunction,x0):
    return Sampler(gridFunction).pointSample( x0 )

def boundarySample(gridFunction, **kwargs):
    return Sampler(gridFunction).boundarySample(**kwargs)
boundarySample.__doc__ = Sampler.boundarySample.__doc__
