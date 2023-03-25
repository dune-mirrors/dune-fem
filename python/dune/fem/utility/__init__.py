from __future__ import print_function
from dune.generator import algorithm, path
import io

# just an idea

def evaluate(expression, coordinate):
    assert isinstance(expression, ufl.Expr)
    return value
def integrate(form,grid=None):
    assert check_form_arity(expression, arguments) == 0
    return value
def assemble(form,space=None):
    arity = check_form_arity(form, arguments)
    assert arity == 1 or arity == 2
    assert space or hasattr(form.testFunction,"space")
    if arity == 1:
        return functional
    else:
        return matrix
def solve(equation,solutionGF):
    arity = check_form_arity(equation, arguments)
    assert arity == 1 or arity == 2
    if arity == 2:
        return info
    else:
        assert solutionGF
        # assert that solutionGF is a coefficient in equation and replace
        # it by trialfunction
        return info

# or all in one?
def evaluate(expression, grid=None, space=None, target=None, coordinate=None):
    if isinstance(expression, ufl.Expr):
        assert coordinate
        return expression(coordinate)
    elif isinstance(expression, ufl.Form):
        if check_form_arity(expression, arguments) == 0:
            assert grid or space
            if not grid: grid = space.gridView
            pass # integrate function
        if check_form_arity(expression, arguments) == 1:
            assert space
            pass # return a df
        elif check_form_arity(expression, arguments) == 2:
            assert space
            pass # return a matrix
    elif isinstance(expression, ufl.Equation):
        assert target
        pass # return solver info

class GridWidth:
    def __init__(self, gridView):
        self._gridView = gridView
        code = """
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/gridpart/common/gridpartadapter.hh>

template <class GridView>
double gridWidth(const GridView& gv)
{
  Dune::Fem::GridPartAdapter< GridView > gp( gv );
  return Dune::Fem::GridWidth::calcGridWidth(gp);
}
"""
        self._gridWidth = algorithm.load("gridWidth", io.StringIO(code), self._gridView )

    def __call__(self):
        """
        Returns:
        --------
            h computed as min(|E|/|e|) forall E in gridView
        """
        return self._gridWidth(self._gridView)

def gridWidth( gridView ):
    """
    Computes the minimal grid width which over all elements by computing min { |E|/|e| }
    for all intersections e of an element.

    Parameters:
        gridView  grid view to compute gridwidth h for. h is comp

    Returns: h computed as min(|E|/|e|) forall E in gridView
    """
    return GridWidth(gridView)()


class Sampler:
    def __init__(self, gridFunction ):
        self.gridFunction = gridFunction
        self.lineSampler = None
        self.pointSampler = None
        self.boundarySampler = None

    def lineSample(self,x0,x1,N):
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
            sample(const GF &gf, DT &start, DT &end, int n)
            {
              Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
              std::vector<DT> coords(n);
              std::vector<typename GF::RangeType> values(n);
              sampler(gf,values);
              sampler.samplePoints(coords);
              return std::make_pair(coords,values);
            }
            #endif

            """
            self.lineSampler = algorithm.load('sample', io.StringIO(_lineSampleCode), self.gridFunction, x0, x1, N)

        p,v = self.lineSampler( self.gridFunction, x0, x1, N )
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


def lineSample(gridFunction,x0,x1,N):
    return Sampler(gridFunction).lineSample(x0, x1, N)

def pointSample(gridFunction,x0):
    return Sampler(gridFunction).pointSample( x0 )

def boundarySample(gridFunction, **kwargs):
    return Sampler(gridFunction).boundarySample(**kwargs)
boundarySample.__doc__ = Sampler.boundarySample.__doc__
