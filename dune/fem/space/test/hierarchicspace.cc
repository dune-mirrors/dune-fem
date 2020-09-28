#include <config.h>

#include <iostream>

using namespace Dune;

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/misc/l2norm.hh>

#include <dune/fem/io/parameter.hh>

#include <dune/fem/test/testgrid.hh>

// polynom approximation order of quadratures,
// at least poolynom order of basis functions
static const int dimRange = 1 ;

//***********************************************************************
/*! check for hierarchical basis
*/
//***********************************************************************

//! the index set we are using
typedef GridSelector::GridType MyGridType;
typedef Dune::Fem::DGAdaptiveLeafGridPart< MyGridType > GridPartType;
typedef Dune::Fem::FunctionSpace < double , double, MyGridType::dimensionworld,
        dimRange> FunctionSpaceType;

typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType,
                                    GridPartType, 0, Dune::Fem::CachingStorage>  DGSpaceType;

typedef Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceType,
                                            GridPartType, 0, Dune::Fem::CachingStorage>  HierarchicLegendreSpaceType;

// ********************************************************************
template <class Space, class SpaceOneOrderMore>
void checkHierarchicSpace ( const Space& space, const SpaceOneOrderMore& spaceOneOrderMore )
{
  typedef typename Space::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > VolumeQuadratureType;

  typedef typename Space::RangeType RangeType;

  std::vector< RangeType > phi;
  std::vector< RangeType > phiMore;

  const int spaceOrder = space.order();
  const int orderMore  = spaceOneOrderMore.order();
  if( spaceOrder+1 != orderMore )
    DUNE_THROW(Dune::InvalidStateException,"order of spaces does not match");

  for( auto it = space.begin(), end = space.end(); it != end; ++it )
  {
    const auto& entity = *it ;

    const auto baseSet     = space.basisFunctionSet( entity );
    const auto baseSetMore = spaceOneOrderMore.basisFunctionSet( entity );

    for( int order = 0; order<2*space.order()+2; ++order )
    {
      VolumeQuadratureType quad( entity, order );

      const int nop = quad.nop();
      for( int qp=0; qp<nop; ++qp )
      {
        const unsigned int size = baseSet.size();
        phi.clear();
        phi.resize( size, RangeType(0) );
        baseSet.evaluateAll( quad[ qp ], phi );

        phiMore.clear();
        phiMore.resize( baseSetMore.size(), RangeType(0) );
        baseSetMore.evaluateAll( quad[ qp ], phiMore );

        if( size >= baseSetMore.size() )
        {
          DUNE_THROW(Dune::InvalidStateException,"Number of basis functions inconsistent");
        }

        bool ok = true;
        for( unsigned int i=0; i<size; ++i )
        {
          if( (phi[ i ] - phiMore[ i ]).two_norm() > 1e-12 )
          {
            std::cerr << std::endl;
            std::cerr << "phi("<< spaceOrder <<")[ " << i << " ] = " << phi[ i ]
                      << "  !=  phi("<<orderMore<< ")[ " << i << " ] = " << phiMore[ i ] << std::endl;
            ok = false ;
          }
        }
        if( ! ok )
          DUNE_THROW(Dune::InvalidStateException,"Basis of V_"<<spaceOrder<<" != V_("<<orderMore<<") / V_(p="<<orderMore<<")");
      }
    }
  }
}

template <class DiscreteSpace,
          int newOrder>
struct ToNewPolorder;

template <class FunctionSpace, class GridPart, int polOrder,
          template <class,class,int,class> class DiscreteSpace,
          int newOrder>
struct ToNewPolorder< DiscreteSpace< FunctionSpace, GridPart, polOrder, Dune::Fem::CachingStorage>, newOrder >
{
  typedef DiscreteSpace< FunctionSpace, GridPart, newOrder, Dune::Fem::CachingStorage> Type;
};


template< class DiscreteSpaceType, bool >
struct CheckSpace
{
  static void check( const DiscreteSpaceType& space )
  {
    static const int polOrder = DiscreteSpaceType :: polynomialOrder ;
    typedef typename ToNewPolorder< DiscreteSpaceType, polOrder+1 >:: Type  NewDiscreteSpaceType;

    NewDiscreteSpaceType spaceMore( space.gridPart() );
    checkHierarchicSpace( space, spaceMore );

    // see orthonormal shapefunction set implementation
    static const int polMax = (DiscreteSpaceType::dimDomain == 2) ? 7 : 3 ;
    CheckSpace< NewDiscreteSpaceType, (polOrder < polMax) > :: check( spaceMore );
  }
};

template< class DiscreteSpaceType >
struct CheckSpace< DiscreteSpaceType, false >
{
  static void check( const DiscreteSpaceType& space )
  {
  }
};

//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main( int argc, char *argv[] )
try {
  Dune::Fem::MPIManager :: initialize( argc, argv );

  const char* paramName = "parameter";
  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std :: endl;
  }
  else
    paramName = argv[1];

  std::string paramFile( paramName );

  // append parameter
  Dune::Fem::Parameter :: append( argc , argv );
  Dune::Fem::Parameter :: append( paramFile );

  if( Dune::Capabilities::hasSingleGeometryType< MyGridType >::v &&
      Dune::GeometryType( Dune::Capabilities::hasSingleGeometryType< MyGridType
        >::topologyId, MyGridType::dimension ).isCube() )
  {
    MyGridType &grid = Dune::Fem::TestGrid::grid();

    GridPartType part ( grid );
    {
      std::cout << "Check DiscontinuousGalerkinSpace... ";
      DGSpaceType space( part );
      CheckSpace< DGSpaceType, Dune::Fem::Capabilities::isHierarchic< DGSpaceType >::v >::check( space );
      std::cout << "successful! " << std::endl;
    }
    {
      std::cout << "Check HierarchicLegendreDiscontinuousGalerkinSpace... ";
      HierarchicLegendreSpaceType space( part );
      CheckSpace< HierarchicLegendreSpaceType,
                  Dune::Fem::Capabilities::isHierarchic< HierarchicLegendreSpaceType >::v >::check( space );
      std::cout << "successful! " << std::endl;
    }
  }
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}

