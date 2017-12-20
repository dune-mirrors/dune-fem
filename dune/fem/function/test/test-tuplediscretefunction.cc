#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

#include <dune/common/hybridutilities.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/tuplediscretefunction.hh>

#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/test/testgrid.hh>

typedef Dune::GridSelector::GridType HGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;

typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, 1 > BaseFunctionSpaceType;

typedef Dune::Fem::DiscontinuousGalerkinSpace< BaseFunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpace1;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< BaseFunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpace2;
typedef Dune::Fem::DiscontinuousGalerkinSpace< BaseFunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpace3;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< BaseFunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpace4;

typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace1 > DiscreteFunction1;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace2 > DiscreteFunction2;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace3 > DiscreteFunction3;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace4 > DiscreteFunction4;

typedef Dune::Fem::TupleDiscreteFunction< DiscreteFunction1, DiscreteFunction2, DiscreteFunction3, DiscreteFunction4 > DiscreteFunctionType;

typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
typedef typename DiscreteFunctionType::Sequence Sequence;

template< class FunctionSpace >
struct MyFunction
{
  typedef FunctionSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = 1.0;
    for( int k = 0; k < FunctionSpaceType::dimDomain; ++k )
      y *= std::sin( M_PI * x[ k ] );
  }

  void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
  {
    for( int j = 0; j < FunctionSpaceType::dimDomain; ++j )
    {
      // jacobian has only one row, calc j-th column
      jacobian[ 0 ][ j ] = M_PI;
      for( int k = 0; k < FunctionSpaceType::dimDomain; ++k )
        jacobian[ 0 ][ j ] *= (j == k ? cos( M_PI*x[ k ] ) : sin( M_PI*x[ k ] ));
    }
    for( int j = 1; j < FunctionSpaceType::dimRange; ++j )
      jacobian[ j ] = jacobian[ 0 ];
  }
};

// main program
int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  try
  {
    auto& grid = Dune::Fem::TestGrid::grid();

    GridPartType gridPart( grid );
    std::cout << "Grid width: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

    DiscreteFunctionSpaceType space( gridPart );

    DiscreteFunctionType df( "ref", space );
    df.clear();

    DiscreteFunctionType df2( "ref2", space );
    df2.clear();

    df2 += df;
    df.axpy( 0.5, df2 );

    std::cout << "dofs = " << df.size() << std::endl;
    std::size_t size( 0 );
    Dune::Hybrid::forEach( Sequence{}, [ & ]( auto i ){ size += std::get< i >( df ).size(); } );
    std::cout << size << std::endl;

    // refine grid
    Dune::Fem::GlobalRefine::apply( grid, 1 );

    df2 += df;
    std::cout << "dofs = " << df.size() << std::endl;
    size = 0;
    Dune::Hybrid::forEach( Sequence{}, [ & ]( auto i ){ size += std::get< i >( df ).size(); } );
    std::cout << size << std::endl;

    std::stringstream stream;
    Dune::Fem::StandardOutStream out( stream );
    df.write( out );
    Dune::Fem::StandardInStream in( stream );
    df.read( in );

    Dune::Fem::L2Norm< GridPartType > l2Norm( gridPart, 5 );

    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FullFunctionSpaceType;
    typedef MyFunction< FullFunctionSpaceType > FullFunction;
    typedef MyFunction< BaseFunctionSpaceType > SubFunction;

    std::cout << "Interpolation tests: " << std::endl;
    for( int i = 0; i < 4; ++i )
    {
      std::cout << "Full function interpolation test:" << std::endl;

      Dune::Fem::interpolate( Dune::Fem::gridFunctionAdapter( FullFunction(), gridPart, 5 ), df );
      const double fullError = l2Norm.distance( df, FullFunction() );
      std::cout<< fullError << std::endl;
      std::cout<< df.name() <<std::endl;

      std::cout << "Checking for subFunctions:" << std::endl;
      Dune::Hybrid::forEach( Sequence{}, [ & ]( auto i ){ std::cout<< std::get< i >( df ).name() << std::endl; } );

      // read access
      Dune::FieldVector< double, 4 > error;
      Dune::Hybrid::forEach( Sequence{}, [ & ]( auto i ){ error[ i ] = l2Norm.distance( std::get< i >( df ), SubFunction() ); } );
      std::cout << error.two_norm() << std::endl;
      assert( std::abs( error.two_norm() - fullError ) < 1e-8 );

      // write access
      Dune::Hybrid::forEach( Sequence{},
        [ & ]( auto i ){ Dune::Fem::interpolate( Dune::Fem::gridFunctionAdapter( SubFunction(), gridPart, 5 ), std::get< i >( df ) ); } );
      std::cout<<std::endl;

      Dune::Fem::GlobalRefine::apply( grid, 1 );
    }

    return 0;
  }
  catch( const Dune::Exception &e )
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
