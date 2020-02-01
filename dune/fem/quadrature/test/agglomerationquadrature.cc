#include <config.h>
#include <iostream>
#include <set>

#include <dune/grid/common/gridinfo.hh>

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/agglomerationquadrature.hh>

#include <dune/fem/quadrature/dunequadratures.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

#include "checkleafcodim1.hh"

using namespace Dune;
using namespace Fem ;

struct DuneQuadratures
{};

template<class GridPartType, bool caching, class QuadratureTraits = int>
class TestCaching
{
private:
  // types of iterators, entities and intersection
  typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
  typedef typename GridPartType :: IndexSetType  IndexSetType ;
  typedef typename IteratorType :: Entity EntityType;
  typedef typename EntityType :: Geometry  Geometry;
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType :: Intersection IntersectionType;

  typedef Dune::Fem::ElementQuadrature<GridPartType, 0>    ElementQuadratureType;

  typedef CachingLumpingQuadrature< GridPartType, 0> LumpingQuadratureType;

  GridPartType & gridPart_;
  int order_;
  const double eps_;
  const bool skipFaces_;

public:

  TestCaching(GridPartType & gridPart, int order, const double eps = 1e-8)
    : gridPart_(gridPart), order_( std::abs(order) ), eps_(eps), skipFaces_( order < 0 )
  { }

  void runTest()
  {
    testElementQuadratures();
  }

  template <class Coordinate>
  double function( const Coordinate& x ) const
  {
    return x * x;
  }

  void testElementQuadratures()
  {
    IteratorType endit = gridPart_.template end<0>();
    for (IteratorType it = gridPart_.template begin<0>(); it != endit; ++it)
    {
      const EntityType & entity = *it;
      const Geometry& geo = entity.geometry();
      const double volume = geo.volume();


      // test integration of 1 over element which should yield the volume.
      double integral = 0 ;
      {
        // test memory management of agglomeration
        {
          ElementQuadratureType quad1( entity, order_, false );
          ElementQuadratureType quad2( entity, order_, false );
        }
        // false here forces the weights and points to be computed
        ElementQuadratureType quad( entity, order_, false );
        const int quadNop = quad.nop();
        // integrate function over element
        for( int qp=0; qp<quadNop; ++qp )
        {
          integral += quad.weight( qp ) * geo.integrationElement( quad.point( qp ) );
        }

        if( std::abs( integral - volume ) > 1e-8 )
        {
          DUNE_THROW(InvalidStateException,"Integral of 1 of element does not yield the volume. Int = " << integral << " vol = " << volume );
        }

        // test integration of function over element and compare with normal
        // quadrature
        integral = 0;
        // integrate function over element
        for( int qp=0; qp<quadNop; ++qp )
        {
          double val = function( geo.global( quad.point( qp ) ) );
          integral += val * quad.weight( qp ) * geo.integrationElement( quad.point( qp ) );
        }
      }

      double checkInt = 0;
      {
        // this only works with regular elements
        assert( ! entity.type().isNone() );
        ElementQuadratureType elemQuad( entity, order_ );
        const int qNop = elemQuad.nop();
        // integrate function over element
        for( int qp=0; qp<qNop; ++qp )
        {
          double val = function( geo.global( elemQuad.point( qp ) ) );
          checkInt += val * elemQuad.weight( qp ) * geo.integrationElement( elemQuad.point( qp ) );
        }
      }

      if( std::abs( integral - checkInt ) > 1e-8 )
      {
        DUNE_THROW(InvalidStateException,"Integral does not equal check! Int = " << integral << " check = " << checkInt );
      }
    }
  }

};

// main program
int main(int argc, char ** argv)
{
  try
  {
    MPIManager :: initialize( argc, argv );
    // *** Initialization
    Parameter::append(argc, argv);
    if (argc==2) Parameter::append(argv[1]);
    else Parameter::append("parameter");

    typedef Dune::GridSelector::GridType GridType;
    std::string filename;
    std::stringstream fileKey ;
    fileKey << "fem.gridfile" << GridType :: dimension ;
    Parameter::get( fileKey.str(), filename);
    Dune::GridPtr< GridType > gridptr( filename );
    GridType &grid = *gridptr;

    int maxlevel = 0; // Parameter::getValue<int>("fem.maxlevel");
    //grid.globalRefine( maxlevel );
    //Dune::gridinfo(grid);

    int quadOrder = 2;//Parameter::getValue<int>("fem.quadorder");

    if ( Parameter::getValue<bool>("fem.skipfaces", false ) )
      quadOrder = -quadOrder ;

    int nonConformOrigin = Parameter::getValue< int > ( "poisson.nonConformOrigin", 0 );
    if ( nonConformOrigin )
    {
      const int refineelement = 1 ;
      std::cout << "Create local refined grid" << std::endl;
      for (int i=0;i<nonConformOrigin;++i)
      {
        if( grid.comm().rank() == 0)
        {
          typedef GridType::Codim<0>::LeafIterator IteratorType;
          IteratorType endit = grid.leafend<0>();
          for(IteratorType it = grid.leafbegin<0>(); it != endit ; ++it)
          {
            const IteratorType :: Entity & entity = *it ;
            const IteratorType :: Entity :: Geometry& geo = entity.geometry();
            if (geo.center().two_norm() < 0.5)
            {
              grid.mark(refineelement, entity );
              std::cout << "mark" << std::endl;
            }
          }
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();
      }
    }

    //typedef HierarchicGridPart< GridType > GridPartType;
    typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
    GridPartType gridPart( grid );

    const double eps = 1e-8;

    for(int l=0; l<=maxlevel; ++l )
    {
      {
        std::cout << "Testing AgglomerationQuadratures: " << std::endl;
        TestCaching<GridPartType,false> testCaching(gridPart, quadOrder, eps);
        testCaching.runTest();
      }
      grid.globalRefine( 1 );
    }
    return 0;
  }
  catch( const Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
