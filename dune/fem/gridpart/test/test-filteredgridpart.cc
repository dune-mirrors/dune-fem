#include<config.h>

// system includes
#include <cassert>
#include <iostream>

// dune-common includes
#include<dune/common/exceptions.hh>
#include <dune/grid/common/genericreferenceelements.hh>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>

// local includes
#include"radialfilter.hh"
#include"../../test/testgrid.hh"


template< class GridPartType >
void testGridPart( const GridPartType & gridPart )
{
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename GridPartType::IndexSetType::IndexType IndexType;
  const IndexSetType & indexSet = gridPart.indexSet();
  IndexType maxIndex = 0;
  int count = 0;
  std::vector< bool > isConsecutiveIndex(GridPartType::dimension+1,true);
  std::vector< std::vector< bool > > consecutiveIndex(GridPartType::dimension+1);
  for (int c=0;c<=GridPartType::dimension;++c)
    consecutiveIndex[c].resize( indexSet.size(c), false );

  typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
  {
    ++count;
    IndexType index = indexSet.index(*it);
    maxIndex = std::max( index, maxIndex);

    if (index >= consecutiveIndex[0].size()) isConsecutiveIndex[0] = false;
    else consecutiveIndex[0][index] = true;
    for (int c=0;c<=GridPartType::dimension;++c)
    {
      int nSubEn = Dune::GenericReferenceElements< typename GridPartType::GridType::ctype, GridPartType::dimension >::
          general( it->type() ).size(c);
      for (int i=0;i<nSubEn;++i)
      {
        IndexType index = indexSet.subIndex(*it,i,c);
        if (index >= consecutiveIndex[c].size()) isConsecutiveIndex[c] = false;
        else consecutiveIndex[c][index] = true;
      }
    }
  }

  std::cout << "entities visited: " << count << std::endl;
  
  std::cout << "entities in index set: " << indexSet.size( 0 ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;

  for (int c=0;c<=GridPartType::dimension;++c)
  {
    std::cout << "index set for codim " << c; 
    if ( !isConsecutiveIndex[c] ) 
    {
      std::cout << " is not consecutive: too large index encountered" << std::endl;
      continue;
    }
    for (unsigned int i=0;i<consecutiveIndex[c].size();++i)
    {
      if (!consecutiveIndex[c][i]) 
      {
        isConsecutiveIndex[c] = false;
        break;
      }
    }
    if ( !isConsecutiveIndex[c] ) 
    {
      std::cout << " is not consecutive: hole encountered" << std::endl;
    }
    else
    {
      std::cout << " is consecutive" << std::endl;
    }
  }
}

template< int codim, class GridPartType >
void testSubEntities( const GridPartType & gridPart )
{
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename GridPartType::IndexSetType::IndexType IndexType;
  const IndexSetType & indexSet = gridPart.indexSet();
  IndexType maxIndex = 0;
  int count = 0;
  typedef typename GridPartType::template Codim< codim >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< codim >();
  for( IteratorType it = gridPart.template begin< codim >(); it != end; ++it )
  {
    ++count;
    IndexType index = indexSet.index(*it);
    maxIndex = std::max( index, maxIndex);
  }

  std::cout << "codim " << codim << " subentities visited: " << count << std::endl;

  std::cout << "entities in index set: " << indexSet.size( codim ) << std::endl;
  std::cout << "maximum value in index set: " << maxIndex << std::endl;
}

typedef Dune::GridSelector::GridType GridType;
typedef Dune::DGAdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > BasicFilterType;
typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, BasicFilterType > FilterType;

int main ( int argc, char ** argv )
{
  Dune::MPIManager :: initialize( argc, argv );
  try
  {
    // create grid
    GridType & grid = Dune::TestGrid::grid();

    // refine grid
    const int step = Dune::TestGrid::refineStepsForHalf();
    grid.globalRefine( 2*step );

    // crete grid part
    HostGridPartType hostGridPart( grid );
    BasicFilterType::GlobalCoordinateType center( 0 );
    BasicFilterType basicFilter( center, .25 );
    FilterType filter( hostGridPart, basicFilter );

    {
      // allow non consecutive index set
      typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, false > GridPartType;
      GridPartType gridPart( hostGridPart, filter );
      // run test
      std::cout << "Non consecutive inex set: test using codim=0 iterator" << std::endl;
      testGridPart( gridPart );
      std::cout << std::endl;
      std::cout << "Non consecutive inex set: test using codim=dimension iterator" << std::endl;
      testSubEntities< GridType::dimension >( gridPart );
      std::cout << std::endl << std::endl;
    }
    {
      // force consecutive index set
      typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, true > GridPartType;
      GridPartType gridPart( hostGridPart, filter );
      // run test
      std::cout << "Consecutive inex set: test using codim=0 iterator" << std::endl;
      testGridPart( gridPart );
      std::cout << std::endl;
      std::cout << "Consecutive inex set: test using codim=dimension iterator" << std::endl;
      testSubEntities< GridType::dimension >( gridPart );
      std::cout << std::endl << std::endl;
    }

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}


