#include <config.h>
#include <iostream>
#include <set>

#include <dune/grid/common/gridinfo.hh>

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/elementquadrature.hh>

#include <dune/fem/quadrature/dunequadratures.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

#include "checkleafcodim1.hh"

using namespace Dune;
using namespace Fem ;

struct DuneQuadratures
{};

template<class GridPartType, int codim, bool caching, class QuadratureTraits = int >
struct QuadratureChooser;

template<class GridPartType, int codim>
struct QuadratureChooser<GridPartType,codim,false,int>
{
  typedef ElementQuadrature< GridPartType, codim > Quadrature;
};
template<class GridPartType, int codim>
struct QuadratureChooser<GridPartType,codim,true,int>
{
  typedef CachingQuadrature< GridPartType, codim > Quadrature;
};

template<class GridPartType, int codim, class QuadratureTraits >
struct QuadratureChooser<GridPartType,codim,false, QuadratureTraits>
{
  typedef ElementQuadrature< GridPartType, codim, Dune::Fem::DuneQuadratureTraits > Quadrature;
};
template<class GridPartType, int codim, class QuadratureTraits >
struct QuadratureChooser<GridPartType,codim,true, QuadratureTraits>
{
  typedef CachingQuadrature< GridPartType, codim, Dune::Fem::DuneQuadratureTraits > Quadrature;
};

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

  typedef typename QuadratureChooser<GridPartType, 0, caching, QuadratureTraits>::Quadrature VolumeQuadratureType;
  typedef typename QuadratureChooser<GridPartType, 0, false,   QuadratureTraits>::Quadrature ElementQuadratureType;
  typedef typename QuadratureChooser<GridPartType, 1, caching>::Quadrature FaceQuadratureType;

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
    if( skipFaces_ )
    {
      std::cout << "Skipping faces due to fem.skipfaces "<<std::endl;
      return ;
    }
    testFaceQuadratures();
  }

  void testElementQuadratures()
  {
    IteratorType endit = gridPart_.template end<0>();
    for (IteratorType it = gridPart_.template begin<0>(); it != endit; ++it)
    {
      const EntityType & entity = *it;
      const Geometry& geo = entity.geometry();

      VolumeQuadratureType cacheQuad( entity, order_ );
      ElementQuadratureType elemQuad( entity, order_ );
      if( cacheQuad.nop() != elemQuad.nop() )
      {
        std::cout << " Error: nops not equal: cachequad = " << cacheQuad.nop()
                  << "   elemQuad = " << elemQuad.nop() << std::endl;
        return ;
      }

      const int quadNop = cacheQuad.nop();
      for (int qp = 0; qp < quadNop; ++qp)
      {
        Dune::FieldVector<typename GridPartType::ctype,GridPartType::dimensionworld> globalElem, globalCache;
        globalCache  = geo.global( cacheQuad.point( qp ) );
        globalElem   = geo.global( elemQuad.point( qp ) );
        if( (globalCache-globalElem).two_norm() > eps_)
        {
          std::cout << " Error: x(cache) = " << globalCache << " != "
                    << globalElem << " = x(elem) " << std::endl;
          break;
        }
      }
    }
  }

  void testFaceQuadratures()
  {
    // just another face check
    CachingQuadratureTest :: checkLeafsCodimOne( gridPart_, order_ );

    const IndexSetType& indexSet = gridPart_.indexSet();

    std::set< int > insideTwists;
    std::set< int > outsideTwists;

    IteratorType endit = gridPart_.template end<0>();
    for (IteratorType it = gridPart_.template begin<0>(); it != endit; ++it)
    {
      const EntityType & entity = *it;
      const IntersectionIteratorType iend = gridPart_.iend(entity);
      for( IntersectionIteratorType iit = gridPart_.ibegin(entity); iit != iend; ++iit )
      {
        const IntersectionType &intersection = *iit;
        if (intersection.neighbor() && intersection.conforming())
        {
          EntityType inside = intersection.inside();
          FaceQuadratureType faceQuadInner(gridPart_, intersection, order_,
              FaceQuadratureType::INSIDE);

          EntityType outside = intersection.outside();
          FaceQuadratureType faceQuadOuter(gridPart_, intersection, order_,
              FaceQuadratureType::OUTSIDE);

          const int faceQuadInner_nop = faceQuadInner.nop();
          for (int qp = 0; qp < faceQuadInner_nop; ++qp)
          {
            typedef typename GridPartType::TwistUtilityType TwistUtilityType;
            Dune::FieldVector<typename GridPartType::ctype,GridPartType::dimensionworld> globalInside, globalOutside;
            globalInside = inside.geometry().global(faceQuadInner.point(qp));
            globalOutside = outside.geometry().global(faceQuadOuter.point(qp));
            if( (globalInside-globalOutside).two_norm() > eps_)
            {
              const int twistInside  = TwistUtilityType::twistInSelf( gridPart_.grid(), intersection );
              const int twistOutside = TwistUtilityType::twistInNeighbor( gridPart_.grid(), intersection );
              insideTwists.insert( twistInside );
              outsideTwists.insert( twistOutside );
              //std::cout << "Inside  twist = " << twistInside << std::endl;
              //std::cout << "Outside twist = " << twistOutside << std::endl;
              std::cout << "On Element " << indexSet.index( entity ) << " with neighbor " << indexSet.index( outside ) << std::endl;
              std::cout << " Error: x(inside) = " << globalInside << " != "
                        << globalOutside << " = x(outside) "
                        << " at intersection.indexInInside() = " << intersection.indexInInside()
                        << ", intersection.indexInOutside() = " << intersection.indexInOutside() << std::endl;
              break;
            }
          }
        }
      }
    }

    typedef typename std::set<int> :: iterator iterator;
    const iterator endin = insideTwists.end();
    for( iterator it = insideTwists.begin(); it != endin; ++it )
    {
      std::cout << "Inside twst: " << (*it) << std::endl;
    }
    const iterator endout = outsideTwists.end();
    for( iterator it = outsideTwists.begin(); it != endout; ++it )
    {
      std::cout << "Outside twst: " << (*it) << std::endl;
    }
  }

  void testLumpingQuadrature()
  {
    IteratorType endit = gridPart_.template end<0>();
    for (IteratorType it = gridPart_.template begin<0>(); it != endit; ++it)
    {
      const EntityType & entity = *it;
      const Geometry& geo = entity.geometry();

      LumpingQuadratureType lumpQuad( entity, order_ );
      const int quadNop = lumpQuad.nop();
      if( quadNop != int(geo.corners()) )
      {
        std::cout << " Error: nops not equal: lumpQuad = " << quadNop
                  << "   geo.corners = " << geo.corners() << std::endl;
        return ;
      }

      for (int qp = 0; qp < quadNop; ++qp)
      {
        Dune::FieldVector<typename GridPartType::ctype,GridPartType::dimensionworld> globalElem, globalCache;
        globalCache  = geo.global( lumpQuad.point( qp ) );
        globalElem   = geo.corner( qp );
        if( (globalCache-globalElem).two_norm() > eps_)
        {
          std::cout << " Error: x(cache) = " << globalCache << " != "
                    << globalElem << " = x(elem) " << std::endl;
          break;
        }
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

    int maxlevel = Parameter::getValue<int>("fem.maxlevel");
    //grid.globalRefine( maxlevel );
    //Dune::gridinfo(grid);

    int quadOrder = Parameter::getValue<int>("fem.quadorder");

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
        std::cout << "Testing default ElementQuadratures: " << std::endl;
        TestCaching<GridPartType,false> testCaching(gridPart, quadOrder, eps);
        testCaching.runTest();
      }
      {
        std::cout << "Testing default CachingQuadratures: " << std::endl;
        TestCaching<GridPartType,true> testCaching(gridPart, quadOrder, eps);
        testCaching.runTest();
      }

      {
        std::cout << "Testing DUNE ElementQuadratures: " << std::endl;
        TestCaching<GridPartType,false, DuneQuadratures> testCaching(gridPart, quadOrder, eps);
        testCaching.runTest();
      }
      {
        std::cout << "Testing DUNE CachingQuadratures: " << std::endl;
        TestCaching<GridPartType,true, DuneQuadratures> testCaching(gridPart, quadOrder, eps);
        testCaching.runTest();
      }
      {
        std::cout << "Testing DUNE CachingLumpiongQuadrature: " << std::endl;
        TestCaching<GridPartType,true, DuneQuadratures> testCaching(gridPart, quadOrder, eps);
        testCaching.testLumpingQuadrature();
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
