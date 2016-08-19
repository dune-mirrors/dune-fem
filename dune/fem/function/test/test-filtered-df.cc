#include <config.h>

#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/filter/threadfilter.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/test/testgrid.hh>

int main( int argc, char** argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  try
  {
    // create grid
    auto& grid = Dune::Fem::TestGrid::grid();

    // refine grid
    const int step = Dune::Fem::TestGrid::refineStepsForHalf();
    grid.globalRefine( 3*step );
    grid.loadBalance();

    // create grid part
    typedef Dune::GridSelector::GridType GridType;
    typedef Dune::Fem::AdaptiveLeafGridPart< GridType > HostGridPartType;
    HostGridPartType hostGridPart( grid );

    // create filters
    typedef std::vector< int > ArrayType;
    ArrayType tags( grid.size(0), 0 );
    auto tagsIt( tags.begin() );
    for( const auto& entity : elements( hostGridPart ) )
    {
      const auto& center( entity.geometry().center() );
      if( center[ 0 ] < 0.5 )
        *tagsIt = 0;
      else
        *tagsIt = 1;
      ++tagsIt;
    }
    typedef Dune::Fem::ThreadFilter< HostGridPartType, ArrayType > FilterType;
    FilterType leftFilter( hostGridPart, tags, 0 );
    FilterType rightFilter( hostGridPart, tags, 1 );

    // create filtered grid parts
    typedef Dune::Fem::FilteredGridPart< HostGridPartType, FilterType, true > FilteredGridPartType;
    FilteredGridPartType leftGridPart( hostGridPart, leftFilter );
    FilteredGridPartType rightGridPart( hostGridPart, rightFilter );

    // check consistency number of entities
    int sizeLeftGridPart( 0 );
    for( auto it = leftGridPart.template begin< 0 >(); it != leftGridPart.template end< 0 >(); ++it )
      ++sizeLeftGridPart;
    int sizeRightGridPart( 0 );
    for( auto it = rightGridPart.template begin< 0 >(); it != rightGridPart.template end< 0 >(); ++it )
      ++sizeRightGridPart;
    if( sizeLeftGridPart+sizeRightGridPart != grid.size( 0 ) )
      DUNE_THROW( Dune::InvalidStateException, "Inconsistent size for grid parts" );

    // check grid parts width
    std::cout << "Host grid part width: " << Dune::Fem::GridWidth::calcGridWidth( hostGridPart )
       << " ( number of entities : " << grid.size( 0 ) << " ) " << std::endl;
    std::cout << "Left grid part width: " << Dune::Fem::GridWidth::calcGridWidth( leftGridPart )
      << " ( number of entities : " << sizeLeftGridPart << " ) " << std::endl;
    std::cout << "Right grid part width: " << Dune::Fem::GridWidth::calcGridWidth( rightGridPart )
      << " ( number of entities : " << sizeRightGridPart << " ) " << std::endl;

    // create dfs and fill dofs
    typedef Dune::Fem::FunctionSpace< double, double, GridType::dimension, 1 > FunctionSpaceType;
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, HostGridPartType, 0 > HostDiscreteFunctionSpaceType;
    HostDiscreteFunctionSpaceType hostSpace( hostGridPart );
    Dune::Fem::AdaptiveDiscreteFunction< HostDiscreteFunctionSpaceType > hostDF ( "host", hostSpace );
    for( const auto& entity : entities( hostDF ) )
    {
      auto localDF = hostDF.localFunction( entity );
      const auto& center( entity.geometry().center() );
      localDF[ 0 ] = ( center[ 0 ] < 0.5 ? center[ 0 ] : center[ 0 ]+2.0 );
    }
    std::cout << "Number of DOFs hostDF : " << hostDF.size() << std::endl;
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, FilteredGridPartType, 0 > FilteredDiscreteFunctionSpaceType;
    FilteredDiscreteFunctionSpaceType leftSpace( leftGridPart );
    Dune::Fem::AdaptiveDiscreteFunction< FilteredDiscreteFunctionSpaceType > leftDF ( "left", leftSpace );
    for( const auto& entity : entities( leftDF ) )
    {
      auto localDF = leftDF.localFunction( entity );
      localDF[ 0 ] = entity.geometry().center()[ 0 ];
    }
    std::cout << "Number of DOFs leftDF : " << leftDF.size() << std::endl;
    FilteredDiscreteFunctionSpaceType rightSpace( rightGridPart );
    Dune::Fem::AdaptiveDiscreteFunction< FilteredDiscreteFunctionSpaceType > rightDF ( "right", rightSpace );
    for( const auto& entity : entities( rightDF ) )
    {
      auto localDF = rightDF.localFunction( entity );
      localDF[ 0 ] = entity.geometry().center()[ 0 ]+2.0;
    }
    std::cout << "Number of DOFs rightDF : " << rightDF.size() << std::endl;

    // check dfs
    for( const auto& entity : entities( leftDF ) )
      if( leftDF.localFunction( entity )[ 0 ] != hostDF.localFunction( entity )[ 0 ] )
        DUNE_THROW( Dune::InvalidStateException, "Inconsistent DOF in leftDF" );
    for( const auto& entity : entities( rightDF ) )
      if( rightDF.localFunction( entity )[ 0 ] != hostDF.localFunction( entity )[ 0 ] )
        DUNE_THROW( Dune::InvalidStateException, "Inconsistent DOF in rightDF" );
  }
  catch( Dune::Exception& e )
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch( ... )
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
