#include <config.h>

#include <iostream>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/test/testgrid.hh>

int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  try
  {
    // create grid
    auto &grid = Dune::Fem::TestGrid::grid();

    // refine grid
    const int step = Dune::Fem::TestGrid::refineStepsForHalf();
    grid.globalRefine( 3*step );
    grid.loadBalance();

    // create grid part
    typedef Dune::GridSelector::GridType GridType;
    typedef Dune::Fem::LeafGridPart< GridType > HostGridPartType;
    HostGridPartType hostGridPart( grid );

    // create filters
    typedef std::vector< int > ArrayType;
    ArrayType tags( grid.size( 0 ), 0 );
    auto tagsIt( tags.begin() );
    for( const auto &entity : elements( hostGridPart ) )
    {
      const auto &center( entity.geometry().center() );
      if( center[ 0 ] < 0.25 )
        *tagsIt = 0;
      else
        *tagsIt = 1;
      ++tagsIt;
    }
    typedef Dune::Fem::DomainFilter< HostGridPartType, ArrayType > FilterType;
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
    typedef Dune::Fem::AdaptiveDiscreteFunction< HostDiscreteFunctionSpaceType > HostDiscreteFunctionType;
    HostDiscreteFunctionType hostDF( "host", hostSpace );
    {
      hostDF.clear();
      Dune::Fem::AddLocalContribution< HostDiscreteFunctionType > localDF( hostDF );
      for( const auto &entity : entities( hostDF ) )
      {
        localDF.bind( entity );
        const auto &center( entity.geometry().center() );
        localDF[ 0 ] = ( center[ 0 ] < 0.25 ? center[ 0 ] : center[ 0 ]+2.0 );
        localDF.unbind();
      }
    }
    std::cout << "Number of DOFs hostDF : " << hostDF.size() << std::endl;
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, FilteredGridPartType, 0 > FilteredDiscreteFunctionSpaceType;
    FilteredDiscreteFunctionSpaceType leftSpace( leftGridPart );
    typedef Dune::Fem::AdaptiveDiscreteFunction< FilteredDiscreteFunctionSpaceType > FilteredDiscreteFunctionType;
    FilteredDiscreteFunctionType leftDF( "left", leftSpace );
    {
      leftDF.clear();
      Dune::Fem::AddLocalContribution< FilteredDiscreteFunctionType > localDF( leftDF );
      for( const auto &entity : entities( leftDF ) )
      {
        localDF.bind( entity );
        localDF[ 0 ] = entity.geometry().center()[ 0 ];
        localDF.unbind();
      }
    }
    std::cout << "Number of DOFs leftDF : " << leftDF.size() << std::endl;
    FilteredDiscreteFunctionSpaceType rightSpace( rightGridPart );
    FilteredDiscreteFunctionType rightDF( "right", rightSpace );
    {
      rightDF.clear();
      Dune::Fem::AddLocalContribution< FilteredDiscreteFunctionType > localDF( rightDF );
      for( const auto &entity : entities( rightDF ) )
      {
        localDF.bind( entity );
        localDF[ 0 ] = entity.geometry().center()[ 0 ]+2.0;
        localDF.unbind();
      }
    }
    std::cout << "Number of DOFs rightDF : " << rightDF.size() << std::endl;

    // check dfs
    Dune::Fem::ConstLocalFunction< HostDiscreteFunctionType > hostLocal( hostDF );
    Dune::Fem::ConstLocalFunction< FilteredDiscreteFunctionType > leftLocal( leftDF );
    Dune::Fem::ConstLocalFunction< FilteredDiscreteFunctionType > rightLocal( rightDF );
    for( const auto &entity : entities( leftDF ) )
    {
      hostLocal.bind( entity );
      leftLocal.bind( entity );
      if( leftLocal[ 0 ] != hostLocal[ 0 ] )
      {
        DUNE_THROW( Dune::InvalidStateException, "Inconsistent DOF in leftDF" );
      }
      hostLocal.unbind();
      leftLocal.unbind();
    }
    for( const auto &entity : entities( rightDF ) )
    {
      hostLocal.bind( entity );
      rightLocal.bind( entity );
      if( rightLocal[ 0 ] != hostLocal[ 0 ] )
        DUNE_THROW( Dune::InvalidStateException, "Inconsistent DOF in rightDF" );
      hostLocal.unbind();
      rightLocal.unbind();
    }

    // create tuple df
    typedef Dune::Fem::TupleDiscreteFunction< FilteredDiscreteFunctionType, FilteredDiscreteFunctionType > TupleDiscreteFunctionType;
    typedef typename TupleDiscreteFunctionType::DiscreteFunctionSpaceType TupleDiscreteSpaceType;
    TupleDiscreteSpaceType tupleSpace(
      std::make_tuple(
        std::make_unique< FilteredDiscreteFunctionSpaceType >( leftGridPart ),
        std::make_unique< FilteredDiscreteFunctionSpaceType >( rightGridPart ) )
      );
    TupleDiscreteFunctionType tupleDF( "tuple", tupleSpace );
    std::cout << "Number of DOFs tupleDF : " << tupleDF.size() << std::endl;
    std::cout << "Number of DOFs tupleDF first component: " << tupleDF.template subDiscreteFunction< 0 >().size() << std::endl;
    std::cout << "Number of DOFs tupleDF second component: " << tupleDF.template subDiscreteFunction< 1 >().size() << std::endl;

    typedef Dune::Fem::TupleDiscreteFunction< FilteredDiscreteFunctionType, HostDiscreteFunctionType > MixedTupleDiscreteFunctionType;
    typedef typename MixedTupleDiscreteFunctionType::DiscreteFunctionSpaceType MixedTupleDiscreteSpaceType;
    MixedTupleDiscreteSpaceType mixedTupleSpace(
      std::make_tuple(
        std::make_unique< FilteredDiscreteFunctionSpaceType >( leftGridPart ),
        std::make_unique< HostDiscreteFunctionSpaceType >( hostGridPart ) )
      );

    MixedTupleDiscreteFunctionType mixedTupleDF( "mixed tuple", mixedTupleSpace );
    std::cout << "Number of DOFs mixedTupleDF : " << mixedTupleDF.size() << std::endl;
    std::cout << "Number of DOFs mixedTupleDF first component: " << mixedTupleDF.template subDiscreteFunction< 0 >().size() << std::endl;
    std::cout << "Number of DOFs mixedTupleDF second component: " << mixedTupleDF.template subDiscreteFunction< 1 >().size() << std::endl;
  }
  catch( Dune::Exception &e )
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
