#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkentityseed.hh>

#include <dune/fem/function/common/common.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpartadapter.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/geogridpart.hh>
#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/lagrange.hh>

#include "failure.hh"
#include "checkseed.hh"
#include "checkgeometry.hh"
#include "checkindexset.hh"
#include "checkintersections.hh"
#include "checkgridpart.hh"
#include <dune/fem/test/testgrid.hh>


template< class GridPart, class LocalFunction >
void testExchangeGeometry ( const GridPart &gridPart, LocalFunction &localFunction )
{
  for( const auto& entity : elements( gridPart ) )
  {
    const auto refElement = Dune::referenceElement< typename GridPart::ctype, GridPart::dimension >( entity.type() );
    localFunction.init( entity.impl().hostEntity() );
    for( int i = 0; i < refElement.size( GridPart::dimension ); ++i )
      for( int k = 0; k < GridPart::dimensionworld; ++k )
        localFunction[ i*GridPart::dimensionworld + k ] = refElement.position( i, GridPart::dimension )[ k ];

    auto xchgEntity = gridPart.exchangeGeometry( entity, localFunction );
    const auto& xchgGeometry = xchgEntity.geometry();
    if( xchgGeometry.type() != entity.type() )
      DUNE_THROW( Dune::InvalidStateException, "exchangeGeometry returns wrong geometry type." );
    if( (xchgGeometry.center() - refElement.position( 0, 0 )).two_norm() > 1e-8 )
    {
      std::cerr << "exchangeGeometry returns wrong center: " << xchgGeometry.center() << std::endl;
      std::cerr << "(real geometry center is: " << entity.geometry().center() << ")." << std::endl;
    }
  }
};

template <class GridPartType>
void testAll( GridPartType& gp )
{
  // test copy constructor
  GridPartType gridPart( gp );

  // test assignment operator and destructor
  {
    GridPartType newGp( gp );
    newGp = gridPart;
  }

  // conversion to and from grid's entity, which is needed for Adaptation
  {
    for( const auto& entity : Dune::elements( gp ) )
    {
      const auto& gridEntity = Dune::Fem::gridEntity( entity );
      const auto& otherEntity = gp.convert( gridEntity );
      otherEntity.type();
    }
  }

  // run tests
  std::cout << "Testing entities" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;

  std::cout << "Testing subentities" << std::endl;
  testSubEntities< GridPartType::dimension >( gridPart );
  std::cout << std::endl;

  std::cout << "GridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  std::cout << "Testing entity seeds" << std::endl;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  std::cout << "Testing geometries" << std::endl;
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing index set" << std::endl;
  Dune::Fem::CheckIndexSet< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing intersections" << std::endl;
  ///Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );

  // check geometries of macro level and leaf level
  //Dune::GeometryChecker<typename GridPartType::GridType> checker;
  //checker.checkGeometry( gridPart );
  Dune::checkGeometryLifetime( gridPart );

  // check entity seeds
  // Dune::checkEntitySeed( gridPart, std::cerr );

  try
  {
    // this test fails if the GridPart is non Cartesian but the grid is
    if( ! Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v )
      checkViewIntersectionIterator( gridPart );

    Dune :: checkIndexSet( gridPart.grid(), gridPart, Dune :: dvverb );
  }
  catch (...)
  {}

  std::cout << std::endl << std::endl;
}

template< bool UseConsecutiveIndexSet, class HostGridPartType, class FilterType >
void testFilteredGridPart( HostGridPartType& hostGridPart, FilterType& filter )
{
  typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, UseConsecutiveIndexSet > GridPartType;
  GridPartType gridPart( hostGridPart, filter );
  testAll( gridPart );
}

template <class HostGridPartType, bool v>
void geomGridPart(HostGridPartType& hostGridPart, std::integral_constant<bool, v> t)
{
  typedef Dune::Fem::FunctionSpace< typename HostGridPartType::ctype, typename HostGridPartType::ctype,
                                    HostGridPartType::dimensionworld, HostGridPartType::dimensionworld > CoordFunctionSpaceType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< CoordFunctionSpaceType, HostGridPartType, 1 > DiscreteCoordFunctionSpaceType;
  DiscreteCoordFunctionSpaceType coordFunctionSpace( hostGridPart );
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteCoordFunctionSpaceType > CoordFunctionType;
  CoordFunctionType coordFunction( "coordinate function", coordFunctionSpace );
  typedef Dune::Fem::Identity< CoordFunctionSpaceType > IdentityType;
  IdentityType identity;
  Dune::Fem::GridFunctionAdapter< IdentityType, HostGridPartType > identitydDF( "identity", identity, hostGridPart );
  Dune::Fem::interpolate( identitydDF, coordFunction );
  typedef Dune::Fem::GeometryGridPart< CoordFunctionType > GridPartType;
  GridPartType gridPart( coordFunction );
  testAll( gridPart );

  if constexpr( v )
  {
    geomGridPart( gridPart, std::false_type() );
  }
}


int main ( int argc, char ** argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  auto& grid = Dune::Fem::TestGrid::grid();

  // refine grid
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();
  grid.globalRefine( 2*step );
  grid.loadBalance();

  typedef  typename Dune::Fem::TestGrid::HGridType GridType;
  // create grid part
  typedef Dune::Fem::LeafGridPart< GridType > HostGridPartType;
  HostGridPartType hostGridPart( grid );

  // AdaptiveLeafGridPart
  {
    std::cout << "*************************************" << std::endl;
    std::cout << "***      AdaptiveLeafridPart      ***"<< std::endl;
    std::cout << "*************************************" << std::endl;

    typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
    typedef typename GridPartType :: GridViewType GridViewType;
    //GridPartType gridPart( grid );
    std::unique_ptr< GridPartType > gpPtr;
    std::unique_ptr< GridViewType > gvPtr;

    gpPtr.reset( new GridPartType( grid ) );
    GridPartType gp2( grid );

    gvPtr.reset( new GridViewType( *gpPtr ) );

    gpPtr.reset( new GridPartType( grid ) );
    gvPtr.reset( new GridViewType( *gpPtr ) );
    gpPtr.reset();
    gvPtr.reset();

    GridPartType& gridPart = gp2;

    testAll( gridPart );
  }

  // GridPartAdapter
  {
    std::cout << "*******************************" << std::endl;
    std::cout << "***    GridPartAdapter      ***"<< std::endl;
    std::cout << "*******************************" << std::endl;

    // create grid part
    typedef Dune::GridSelector::GridType GridType;
    typedef Dune::Fem::GridPartAdapter< typename GridType::LeafGridView >  GridPartType;

    auto gv = grid.leafGridView();
    // GridPartAdapter only stores a pointer
    GridPartType gridPart( gv );
    testAll( gridPart );
  }


  // IdGridPart
  {
    std::cout << "*******************************" << std::endl;
    std::cout << "***      IdGridPart         ***"<< std::endl;
    std::cout << "*******************************" << std::endl;

    typedef Dune::Fem::IdGridPart< HostGridPartType > GridPartType;
    GridPartType gridPart( hostGridPart );

    testAll( gridPart );
  }

  // GeoGridPart
  {
    std::cout << "*******************************" << std::endl;
    std::cout << "***      GeoGridPart        ***"<< std::endl;
    std::cout << "*******************************" << std::endl;

    typedef Dune::Fem::FunctionSpace< GridType::ctype, GridType::ctype, GridType::dimensionworld, GridType::dimensionworld > CoordFunctionSpaceType;
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< CoordFunctionSpaceType, HostGridPartType, 1 > DiscreteCoordFunctionSpaceType;
    DiscreteCoordFunctionSpaceType coordFunctionSpace( hostGridPart );
    typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteCoordFunctionSpaceType > CoordFunctionType;
    CoordFunctionType coordFunction( "coordinate function", coordFunctionSpace );
    typedef Dune::Fem::Identity< CoordFunctionSpaceType > IdentityType;
    IdentityType identity;
    Dune::Fem::GridFunctionAdapter< IdentityType, HostGridPartType > identitydDF( "identity", identity, hostGridPart );
    Dune::Fem::interpolate( identitydDF, coordFunction );
    typedef Dune::Fem::GeoGridPart< CoordFunctionType > GridPartType;
    GridPartType gridPart( coordFunction );

    testAll( gridPart );

    std::cout << "Testing exchange geometry" << std::endl;
    Dune::Fem::TemporaryLocalFunction< DiscreteCoordFunctionSpaceType > tlf( coordFunctionSpace );
    testExchangeGeometry( gridPart, tlf );
    std::cout << std::endl;
  }

  // GeometryGridPart
  {
    std::cout << "************************************" << std::endl;
    std::cout << "***      GeometryGridPart        ***"<< std::endl;
    std::cout << "************************************" << std::endl;

    geomGridPart( hostGridPart, std::true_type() );
  }

  // FilteredGridPart
  {
    std::cout << "******************************" << std::endl;
    std::cout << "***    FilteredGridPart    ***"<< std::endl;
    std::cout << "******************************" << std::endl;

    // create radial filter
    typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > RadialFilterType;
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, RadialFilterType > WrapperRadialFilterType;
    WrapperRadialFilterType wrappedRadialFilter( hostGridPart, RadialFilterType::GlobalCoordinateType( 0 ), 0.25 );

    // test FilteredGridPart with radial filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );

    // test FilteredGridPart with radial filter and forcing consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: force consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< true, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );
  }


  {
    std::cout << "************************************" << std::endl;
    std::cout << "***    FilteredGridPart (tags)   ***"<< std::endl;
    std::cout << "************************************" << std::endl;

    // create domain filter
    typedef std::vector< int > DomainArrayType;
    DomainArrayType tags( grid.size(0), 0 );
    for( std::size_t i = ( tags.size()/2 ); i < tags.size(); ++i )
      tags[ i ] = 1;
    typedef Dune::Fem::DomainFilter< HostGridPartType, DomainArrayType > DomainFilterType;
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, DomainFilterType > WrapperDomainFilterType;
    WrapperDomainFilterType wrapperDomainFilter( hostGridPart, hostGridPart, tags, 1 );

    // test FilteredGridPart with domain filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with domain filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, WrapperDomainFilterType >( hostGridPart, wrapperDomainFilter );
  }

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
