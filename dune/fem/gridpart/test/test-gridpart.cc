#include<config.h>

//- system includes
#include <cassert>

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
  IndexSetType indexSet = gridPart.indexSet();
  std::cout << "entities in index set: " << indexSet.size( 0 ) << std::endl;
}


typedef Dune::GridSelector::GridType GridType;
typedef Dune::DGAdaptiveLeafGridPart< GridType > ContainedGridPartType;
#if defined  USE_FILTEREDGRIDPART
typedef Dune::Fem::RadialFilter< ContainedGridPartType > FilterType;
typedef Dune::Fem::FilteredGridPart< FilterType, true > GridPartType;
#else
typedef ContainedGridPartType GridPartType;
#endif

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
    GridPartType gridPart( grid );

    // run test
    testGridPart( gridPart );

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}


