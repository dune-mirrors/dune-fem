#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
// do nothing in this test if experimental_grid_extension is not activated
int main () { return 0; }
#else

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
#include <dune/fem/test/testgrid.hh>

template< class GridPartType >
void testGridPart( const GridPartType & gridPart )
{
  const auto& indexSet = gridPart.indexSet();
  typedef typename GridPartType::IndexSetType::IndexType IndexType;
  IndexType maxIndex = 0;
  int count = 0;
  std::vector< bool > isConsecutiveIndex(GridPartType::dimension+1,true);
  std::vector< std::vector< bool > > consecutiveIndex(GridPartType::dimension+1);
  for(int c=0;c<=GridPartType::dimension;++c)
    consecutiveIndex[c].resize( indexSet.size(c), false );

  for( const auto& element : elements( gridPart ) )
  {
    ++count;
    auto index = indexSet.index(element);
    maxIndex = std::max( index, maxIndex);

    if(index >= IndexType(consecutiveIndex[0].size()))
      isConsecutiveIndex[0] = false;
    else
      consecutiveIndex[0][index] = true;
    for(int c=0;c<=GridPartType::dimension;++c)
    {
      int nSubEn =
        Dune::referenceElement< typename GridPartType::GridType::ctype, GridPartType::dimension >( element.type() ).size(c);
      for(int i=0;i<nSubEn;++i)
      {
        auto index = indexSet.subIndex(element,i,c);
        if(index >= IndexType(consecutiveIndex[c].size()))
          isConsecutiveIndex[c] = false;
        else
          consecutiveIndex[c][index] = true;
      }
    }
  }

  std::cout << "entities visited: " << count << std::endl;
  std::cout << "entities in index set: " << indexSet.size( 0 ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;

  for(int c=0;c<=GridPartType::dimension;++c)
  {
    std::cout << "index set for codim " << c;
    if( !isConsecutiveIndex[c] )
    {
      std::cout << " is not consecutive: too large index encountered" << std::endl;
      continue;
    }
    for(unsigned int i=0;i<consecutiveIndex[c].size();++i)
    {
      if(!consecutiveIndex[c][i])
      {
        isConsecutiveIndex[c] = false;
        break;
      }
    }
    if( !isConsecutiveIndex[c] )
      std::cout << " is not consecutive: hole encountered" << std::endl;
    else
      std::cout << " is consecutive" << std::endl;
  }
}

template< int codim, class GridPartType >
void testSubEntities( const GridPartType & gridPart )
{
  const auto& indexSet = gridPart.indexSet();
  typename GridPartType::IndexSetType::IndexType maxIndex = 0;
  int count = 0;
  for( auto it = gridPart.template begin< codim >(); it != gridPart.template end< codim >(); ++it )
  {
    ++count;
    auto index = indexSet.index(*it);
    maxIndex = std::max( index, maxIndex);
  }

  std::cout << "codim " << codim << " subentities visited: " << count << std::endl;
  std::cout << "entities in index set: " << indexSet.size( codim ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;
}

template< class GridPartType >
void testIntersectionIterator( const GridPartType & gridPart )
{
  std::vector<int> index( gridPart.indexSet().size(0), 0 );

  for( const auto& element : elements( gridPart ) )
    index[ gridPart.indexSet().index( element ) ] = 1;

  for( const auto& element : elements( gridPart ) )
    for( const auto& intersection : intersections( gridPart, element ) )
      if( intersection.neighbor() )
      {
        const auto& neighbor = intersection.outside();
        const auto& nbIndex = gridPart.indexSet().index( neighbor );
        if( nbIndex >= index.size() )
        {
          std::cout << "An index on neighbor is too large" << std::endl;
          continue;
        }
        if( index[ nbIndex ] == 0 )
        {
          std::cout << "A neighbor is not part of the gridPart" << std::endl;
          continue;
        }
      }
}


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
#endif // #if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
