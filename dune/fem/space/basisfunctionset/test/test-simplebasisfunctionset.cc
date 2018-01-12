#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/function/common/functionset.hh>
#include <dune/fem/function/localfunction/localfunctionsetadapter.hh>
#include <dune/fem/space/fourier/functionset.hh>

#include <dune/fem/space/basisfunctionset/simple.hh>

#include "checkbasisfunctionset.hh"
#include <dune/fem/test/testgrid.hh>


template< class GridPartType, int polorder >
void traverse ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, 1 > ScalarFunctionSpaceType;

  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  auto iterator = gridPart.template begin< 0 >();
  iterator++;
  const EntityType &entity = *iterator;

  // create quadrature
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  QuadratureType quadrature( entity, polorder );

  // type of error
  typedef Dune::FieldVector< double, 7 > ErrorType;

  double eps = 1e-8;

  ErrorType error( 0 );

  typedef Dune::Fem::FourierFunctionSet< ScalarFunctionSpaceType, polorder > FunctionSetType;
  typedef Dune::Fem::FunctionSetProxy< FunctionSetType > FunctionSetProxyType;
  typedef Dune::Fem::LocalFunctionSetAdapter< EntityType, FunctionSetProxyType > LocalFunctionSetType;

  FunctionSetType functionSet( polorder );

  Dune::Fem::SimpleBasisFunctionSet< LocalFunctionSetType >  basisSet( LocalFunctionSetType( entity, &functionSet ) );
  error = Dune::Fem::checkQuadratureConsistency( basisSet, quadrature, true );
  if( error.two_norm() > eps )
  {
    std::cerr<<"Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy ): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, "SimpleBasisFunctionSet test failed." );
  }
}


int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( argc >= 2 ? argv[ 1 ] : "parameter" );

  typedef Dune::GridSelector::GridType GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();

  grid.globalRefine( 1 );

  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  if( gridPart.begin< 0 >() == gridPart.end< 0 >() )
    return 1;

  traverse< GridPartType, 1 >( gridPart );
  traverse< GridPartType, 2 >( gridPart );
  traverse< GridPartType, 3 >( gridPart );
  traverse< GridPartType, 4 >( gridPart );

  return 0;
}
