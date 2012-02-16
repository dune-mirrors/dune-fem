#include <config.h>

#include <cassert>
#include <iostream>

#include <dune/common/exceptions.hh>

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/referenceelements.hh>
#else
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/hierarchicgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

#include "./failure.hh"
#include "./checkseed.hh"
#include "../../test/testgrid.hh"

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

template< class GridPartType >
void testIntersectionIterator( const GridPartType & gridPart )
{
  std::vector<int> index( gridPart.indexSet().size(0), 0 );
  typedef typename GridPartType::template Codim< 0 >::template Partition<Dune::All_Partition>::IteratorType IteratorType;
  const IteratorType end = gridPart.template end< 0,Dune::All_Partition >();
  for( IteratorType it = gridPart.template begin< 0,Dune::All_Partition >(); it != end; ++it )
    index[ gridPart.indexSet().index( * it ) ] = 1;
  for( IteratorType it = gridPart.template begin< 0,Dune::All_Partition >(); it != end; ++it )
  {
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    const IntersectionIteratorType iend = gridPart.iend( *it );
    for ( IntersectionIteratorType inter = gridPart.ibegin( *it );
          inter != iend; ++inter )
    {
      if (inter->neighbor())
      {
        typename GridPartType::IndexSetType::IndexType nbIndex = gridPart.indexSet().index( *(inter->outside()) );
        if ( nbIndex >= index.size() )
        {
          std::cout << "An index on neighbor is too large" << std::endl;
          continue;
        }
        if ( index[ nbIndex ] == 0 )
        {
          std::cout << "A neighbor is not part of the gridPart" << std::endl;
          continue;
        }
      }
    }
  }
}



typedef Dune::GridSelector::GridType GridType;
typedef Dune::DGAdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::IdGridPart< HostGridPartType > GridPartType;

int main ( int argc, char ** argv )
try
{
  Dune::MPIManager::initialize( argc, argv );

  // create grid
  GridType &grid = Dune::TestGrid::grid();

  // refine grid
  const int step = Dune::TestGrid::refineStepsForHalf();
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
  std::cout << "gridWidth: " << Dune::GridWidth::calcGridWidth( gridPart ) << std::endl;

  std::cout << std::endl << std::endl;

  // check entity seed
  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  Dune::Fem::Check< GridPartType, FailureHandlerType >::checkEntitySeed( gridPart, failureHandler );

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
