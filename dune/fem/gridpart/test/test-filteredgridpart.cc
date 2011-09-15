#include<config.h>

// system includes
#include <cassert>
#include <iostream>

// dune-common includes
#include<dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>

// local includes
#include"radialfilter.hh"
#include"../../test/testgrid.hh"


template< class GridPartType >
void testGridPart( const GridPartType & gridPart )
{
  int count = 0;
  typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
    ++count;

  std::cout << "entities visited: " << count << std::endl;
  
  typedef typename GridPartType::IndexSetType IndexSetType;
  const IndexSetType & indexSet = gridPart.indexSet();
  std::cout << "entities in index set: " << indexSet.size( 0 ) << std::endl;
}

template< int codim, class GridPartType >
void testSubEntities( const GridPartType & gridPart )
{
  int count = 0;
  typedef typename GridPartType::template Codim< codim >::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< codim >();
  for( IteratorType it = gridPart.template begin< codim >(); it != end; ++it )
    ++count;

  std::cout << "codim " << codim << " subentities visited: " << count << std::endl;

  typedef typename GridPartType::IndexSetType IndexSetType;
  const IndexSetType & indexSet = gridPart.indexSet();
  std::cout << "entities in index set: " << indexSet.size( codim ) << std::endl;

}


typedef Dune::GridSelector::GridType GridType;
typedef Dune::DGAdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > BasicFilterType;
typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, BasicFilterType > FilterType;
typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, true > GridPartType;

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
    GridPartType gridPart( hostGridPart, filter );

    // run test
    testGridPart( gridPart );
    testSubEntities< GridType::dimension >( gridPart );

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}


