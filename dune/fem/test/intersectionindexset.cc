#undef NDEBUG

// without this define this code won't compile for Yasp
#define ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/hierarchicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

#if defined  USE_FILTEREDGRID 
#include <dune/fem/gridpart/filteredgrid.hh>
#endif

#include "testgrid.hh"
#include "dfspace.hh"
#include "exactsolution.hh"

using namespace Dune;

typedef GridSelector::GridType MyGridType;
// typedef HierarchicGridPart< MyGridType >  ContainedGridPartType;
typedef IntersectionAdaptiveLeafGridPart< MyGridType > ContainedGridPartType;

// use filtered grid for testing 
#if defined  USE_FILTEREDGRID 
  typedef RadialFilter< ContainedGridPartType > FilterType;
  typedef FilteredGridPart<ContainedGridPartType, FilterType, true > GridPartType;
#else
  typedef ContainedGridPartType GridPartType;
#endif

template <class GridPartType> 
void checkIntersectionIndexSet( const GridPartType& gridPart ) 
{
  enum { dim = GridPartType :: GridType :: dimension };
  typedef typename GridPartType :: IndexSetType IndexSetType ;
  const IndexSetType& indexSet = gridPart.indexSet();

  typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType ;
  typedef typename IteratorType :: Entity  EntityType ;
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef typename GridPartType :: IntersectionType         IntersectionType;
  
  int count = 0 ;
  const IteratorType endit = gridPart.template end<0> ();
  for( IteratorType it = gridPart.template begin<0> ();
       it != endit; ++it ) 
  {
    const EntityType& entity = *it;
    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for( IntersectionIteratorType iit = gridPart.ibegin( entity );
         iit != endiit ; ++ iit ) 
    {
      const IntersectionType&  intersection = *iit ; 

      const int index = indexSet.index( *iit ); 
      std::cout << index << " index of intersection " << std::endl;
      ++ count ;
      if( intersection.boundary() ) 
        ++ count ;
    }
  }

  const size_t noIntersections = count / 2;
  std::cout << indexSet.size( dim + 1 ) << " size " << std::endl;
  if( noIntersections != indexSet.size( dim + 1 ) ) 
  {
    std::cerr << "ERROR: size of intersection set wrong! " << std::endl;
  }
}

// main program 
int main(int argc, char ** argv) 
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );

    GridPartType gridPart( grid );
    // add check for grid width 
    std::cout << "Grid width: " 
      << GridWidth :: calcGridWidth( gridPart ) << std::endl; 

    checkIntersectionIndexSet( gridPart );
    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
