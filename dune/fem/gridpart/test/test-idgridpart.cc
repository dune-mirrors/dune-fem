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

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

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
  typename GridPartType::IndexSetType::IndexType maxIndex = 0;
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

    if(static_cast<std::size_t>(index) >= consecutiveIndex[0].size())
      isConsecutiveIndex[0] = false;
    else
      consecutiveIndex[0][index] = true;
    for(int c=0;c<=GridPartType::dimension;++c)
    {
      const int nSubEn =
        Dune::ReferenceElements< typename GridPartType::GridType::ctype, GridPartType::dimension >::general( element.type() ).size(c);
      for(int i=0;i<nSubEn;++i)
      {
        const auto index = indexSet.subIndex(element,i,c);
        if(static_cast<std::size_t>(index) >= consecutiveIndex[c].size())
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

  for( const auto& element : elements( gridPart )  )
    for( const auto& intersection : intersections( gridPart, element ) )
      if( intersection.neighbor() )
      {
        const auto neighbor = intersection.outside();
        const auto nbIndex = gridPart.indexSet().index( neighbor );
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



typedef Dune::GridSelector::GridType GridType;
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::IdGridPart< HostGridPartType > GridPartType;

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
  HostGridPartType hostGridPart( grid );
  GridPartType gridPart( hostGridPart );

  // run test
  std::cout << "test using codim=0 iterator" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;
  std::cout << "test using codim=dimension iterator" << std::endl;
  testSubEntities< GridType::dimension >( gridPart );

  std::cout << std::endl;
  std::cout << "gridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  std::cout << std::endl << std::endl;

  // check entity seed
  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  Dune::Fem::CheckIndexSet< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
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
