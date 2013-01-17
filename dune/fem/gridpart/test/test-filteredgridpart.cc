#include<config.h>

//- system includes
#include <cassert>
#include <iostream>

//- dune-common includes
#include<dune/common/exceptions.hh>

//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>

//- dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/misc/gridwidth.hh>

//- dune-fem test includes
#include <dune/fem/gridpart/test/failure.hh>
#include <dune/fem/gridpart/test/checkseed.hh>
#include <dune/fem/gridpart/test/checkgeometry.hh>
#include <dune/fem/gridpart/test/checkintersections.hh>
#include <dune/fem/test/testgrid.hh>


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
      int nSubEn = Dune::ReferenceElements< typename GridPartType::GridType::ctype, GridPartType::dimension >::
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
typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType > HostGridPartType;
typedef Dune::Fem::RadialFilter< GridType::ctype, GridType::dimensionworld > BasicFilterType;
typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, BasicFilterType > FilterType;

int main ( int argc, char ** argv )
{
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try
  {
    // create grid
    GridType & grid = Dune::Fem::TestGrid::grid();

    // refine grid
    const int step = Dune::Fem::TestGrid::refineStepsForHalf();
    grid.globalRefine( 2*step );
    grid.loadBalance();

    // create grid part
    HostGridPartType hostGridPart( grid );
    BasicFilterType::GlobalCoordinateType center( 0 );
    BasicFilterType basicFilter( center, .25 );
    FilterType filter( hostGridPart, basicFilter );

    {
      // allow non consecutive index set
      typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, false > GridPartType;
      GridPartType gridPart( hostGridPart, filter );
      // run test
      std::cout << "Non consecutive index set: test using codim=0 iterator" << std::endl;
      testGridPart( gridPart );
      std::cout << std::endl;
      std::cout << "Non consecutive index set: test using codim=dimension iterator" << std::endl;
      testSubEntities< GridType::dimension >( gridPart );

      /* ----------- added to test the intersection iterator ---------------*/
      std::cout << std::endl;
      std::cout << "gridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;
      /* --------------------------*/

      /* ----------- new test suite methods ---------------*/
      std::cout << std::endl;
      typedef Dune::DefaultFailureHandler FailureHandlerType;
      FailureHandlerType failureHandler;
      std::cout << "Testing entity seeds" << std::endl;
      Dune::Fem::CheckEntitySeed< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      std::cout << "Testing geometies" << std::endl;
      Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      std::cout << "Testing intersections" << std::endl;
      Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      /* --------------------------*/

      std::cout << std::endl << std::endl;
    }
    {
      // force consecutive index set
      typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, true > GridPartType;
      GridPartType gridPart( hostGridPart, filter );
      // run test
      std::cout << "Consecutive index set: test using codim=0 iterator" << std::endl;
      testGridPart( gridPart );
      std::cout << std::endl;
      std::cout << "Consecutive index set: test using codim=dimension iterator" << std::endl;
      testSubEntities< GridType::dimension >( gridPart );

      /* ----------- added to test the intersection iterator ---------------*/
      std::cout << std::endl;
      std::cout << "gridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;
      /* --------------------------*/

      /* ----------- new test suite methods ---------------*/
      std::cout << std::endl;
      typedef Dune::DefaultFailureHandler FailureHandlerType;
      FailureHandlerType failureHandler;
      std::cout << "Testing entity seeds" << std::endl;
      Dune::Fem::CheckEntitySeed< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      std::cout << "Testing geometies" << std::endl;
      Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      std::cout << "Testing intersections" << std::endl;
      Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
      /* --------------------------*/

      std::cout << "Testing intersection iterator" << std::endl;
      testIntersectionIterator( gridPart );

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
