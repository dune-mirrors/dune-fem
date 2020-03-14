#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/function/common/common.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/geogridpart.hh>
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

  // create grid part
  typedef Dune::GridSelector::GridType GridType;
  typedef Dune::Fem::LeafGridPart< GridType > HostGridPartType;
  HostGridPartType hostGridPart( grid );
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

  // run tests
  std::cout << "Testing entities" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;

  std::cout << "Testing subentities" << std::endl;
  testSubEntities< GridType::dimension >( gridPart );
  std::cout << std::endl;

  std::cout << "GridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  std::cout << "Testing exchange geometry" << std::endl;
  Dune::Fem::TemporaryLocalFunction< DiscreteCoordFunctionSpaceType > tlf( coordFunctionSpace );
  testExchangeGeometry( gridPart, tlf );
  std::cout << std::endl;

  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  std::cout << "Testing entity seeds" << std::endl;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  std::cout << "Testing geometries" << std::endl;
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing index set" << std::endl;
  Dune::Fem::CheckIndexSet< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing intersections" << std::endl;
  Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );

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
