#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/threadfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/misc/gridwidth.hh>

#include "failure.hh"
#include "checkseed.hh"
#include "checkgeometry.hh"
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
      const auto nSubEn =
        Dune::ReferenceElements< typename GridPartType::GridType::ctype, GridPartType::dimension >::general( element.type() ).size(c);
      for(int i=0;i<nSubEn;++i)
      {
        const auto& index = indexSet.subIndex(element,i,c);
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
    const auto& entity = *it;
    ++count;
    const auto index = indexSet.index(entity);
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

  for( const auto& element : elements( gridPart ))
    for( const auto& intersection : intersections( gridPart, element ) )
      if(intersection.neighbor())
      {
        const auto nbIndex = gridPart.indexSet().index( intersection.outside() );
        if( static_cast<std::size_t>(nbIndex) >= index.size() )
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

template< bool UseConsecutiveIndexSet, class HostGridPartType, class FilterType >
void testFilteredGridPart( HostGridPartType& hostGridPart, FilterType& filter )
{
  typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, UseConsecutiveIndexSet > GridPartType;
  GridPartType gridPart( hostGridPart, filter );

  std::cout << "Testing entities" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;

  std::cout << "Testing subentities" << std::endl;
  testSubEntities< HostGridPartType::GridType::dimension >( gridPart );
  std::cout << std::endl;

  std::cout << "GridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  std::cout << "Testing entity seeds" << std::endl;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  std::cout << "Testing geometies" << std::endl;
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing intersections" << std::endl;
  Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
}

int main( int argc, char ** argv )
{
  Dune::Fem::MPIManager :: initialize( argc, argv );

  try
  {
    // create grid
    auto& grid = Dune::Fem::TestGrid::grid();

    // refine grid
    const int step = Dune::Fem::TestGrid::refineStepsForHalf();
    grid.globalRefine( 2*step );
    grid.loadBalance();

    // create grid part
    typedef Dune::GridSelector::GridType GridType;
    typedef Dune::Fem::AdaptiveLeafGridPart< GridType > HostGridPartType;
    HostGridPartType hostGridPart( grid );

    // create radial filter
    typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > RadialFilterType;
    RadialFilterType::GlobalCoordinateType center( 0 );
    RadialFilterType radialFilter( center, .25 );
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, RadialFilterType > WrapperRadialFilterType;
    WrapperRadialFilterType wrappedRadialFilter( hostGridPart, radialFilter );

    // test FilteredGridPart with radial filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );

    // test FilteredGridPart with radial filter and forcing consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with radial filter: force consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< true, HostGridPartType, WrapperRadialFilterType >( hostGridPart, wrappedRadialFilter );

    // create thread filter
    typedef std::vector< int > ThreadArrayType;
    ThreadArrayType tags( grid.size(0), 0 );
    for( std::size_t i = ( tags.size()/2 ); i < tags.size(); ++i )
      tags[ i ] = 1;
    typedef Dune::Fem::ThreadFilter< HostGridPartType, ThreadArrayType > ThreadFilterType;
    ThreadFilterType threadFilter( hostGridPart, tags, 1 );

    // test FilteredGridPart with thread filter and allowing non consecutive index set
    std::cout << std::endl << "Testing FilteredGridPart with thread filter: allow non consecutive index set" << std::endl << std::endl;
    testFilteredGridPart< false, HostGridPartType, ThreadFilterType >( hostGridPart, threadFilter );
  }
  catch(Dune::Exception &e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
