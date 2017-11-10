#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/test/testgrid.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/tuplediscretefunction.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/operator/common/tuple.hh>

#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2norm.hh>

#include "simpleoperators.hh"


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

typedef Dune::Fem::TupleDiscreteFunction< DiscreteFunction1, DiscreteFunction2, DiscreteFunction3, DiscreteFunction4 >
DiscreteFunctionType;

typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

typedef Dune::Fem::TupleOperator<
  SimpleMassOperator< DiscreteFunction1, DiscreteFunction1 >,
  NullOp< DiscreteFunction2, DiscreteFunction1 >,
  NullOp< DiscreteFunction3, DiscreteFunction1 >,
  NullOp< DiscreteFunction4, DiscreteFunction1 > > RowOperator1;

typedef Dune::Fem::TupleOperator<
  NullOp< DiscreteFunction1, DiscreteFunction2 >,
  SimpleMassOperator< DiscreteFunction2, DiscreteFunction2 >,
  NullOp< DiscreteFunction3, DiscreteFunction2 >,
  SimpleLaplaceOperator< DiscreteFunction4, DiscreteFunction2 > > RowOperator2;

typedef Dune::Fem::TupleOperator<
  NullOp< DiscreteFunction1, DiscreteFunction3 >,
  NullOp< DiscreteFunction2, DiscreteFunction3 >,
  NullOp< DiscreteFunction3, DiscreteFunction3 >,
  SimpleMassOperator< DiscreteFunction4, DiscreteFunction3 > > RowOperator3;

typedef Dune::Fem::TupleOperator<
  NullOp< DiscreteFunction1, DiscreteFunction4 >,
  SimpleLaplaceOperator< DiscreteFunction2, DiscreteFunction4 >,
  SimpleMassOperator< DiscreteFunction3, DiscreteFunction4 >,
  NullOp< DiscreteFunction4, DiscreteFunction4 > > RowOperator4;

typedef Dune::Fem::RowTupleOperator< RowOperator1, RowOperator2, RowOperator3, RowOperator4 > OperatorType;


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
    HGridType &grid = Dune::Fem::TestGrid::grid();

    GridPartType gridPart( grid );
    // add check for grid width
    std::cout << "Grid width: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;

    DiscreteFunctionSpaceType space( gridPart );

    DiscreteFunctionType arg( "ref", space );

    DiscreteFunctionType dest( "ref2", space );

    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FullFunctionSpaceType;
    typedef MyFunction< FullFunctionSpaceType > FullFunction;
    Dune::Fem::interpolate( Dune::Fem::gridFunctionAdapter( FullFunction(), gridPart, 5 ), arg );

    std::string name = "bla";
    int second = -42;

    OperatorType mass(
      std::make_tuple(
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ) ),
      std::make_tuple(
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ) ),
      std::make_tuple(
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ) ),
      std::make_tuple(
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ),
        std::tie( name, second ) )
      );

    std::cout<<"Interpolation tests: "<<std::endl;
    for( int i = 0; i < 4; ++i )
    {
      mass( arg, dest );

      std::get< 0 >( mass ) ( arg, std::get< 0 >( dest ) );
      std::get< 2 >( std::get< 0 >( mass ) )( std::get< 2 >( arg ), std::get< 0 >( dest ) );

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
