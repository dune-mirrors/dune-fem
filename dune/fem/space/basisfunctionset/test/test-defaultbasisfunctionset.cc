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

#include <dune/fem/space/shapefunctionset/lagrange.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/simple.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include "checkbasisfunctionset.hh"
#include <dune/fem/test/testgrid.hh>


template< class GridPartType, int polorder >
void traverse ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimRange = DIMRANGE;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, 1 > ScalarFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, dimRange > FunctionSpaceType;

  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  auto iterator = gridPart.template begin< 0 >();
  iterator++;
  const EntityType &entity = *iterator;

  // create quadrature
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  QuadratureType quadrature( entity, polorder );

  // needs a geometry type to construct
  typedef Dune::Fem::LagrangeShapeFunctionSet< ScalarFunctionSpaceType, polorder > ScalarLagrangeShapeFunctionSetType;

  // needs an order to construct
  typedef Dune::Fem::LegendreShapeFunctionSet< ScalarFunctionSpaceType > ScalarLegendreShapeFunctionSetType;

  // needs a geometry type to construct
  typedef Dune::Fem::OrthonormalShapeFunctionSet< ScalarFunctionSpaceType, polorder > ScalarOrthonormalShapeFunctionSetType;

  // type of error
  typedef Dune::FieldVector< double, 5 > ErrorType;

  // prepare shapefunctions
  ScalarLagrangeShapeFunctionSetType scalarLagrangeShapeFunctionSet( entity.type() );
  ScalarLegendreShapeFunctionSetType scalarLegendreShapeFunctionSet( polorder );
  ScalarOrthonormalShapeFunctionSetType scalarOrthonormalShapeFunctionSet( entity.type() );

  double eps = 1e-8;

  ErrorType error( 0 );
  // default basis function set
  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLagrangeShapeFunctionSetType >
  basisSet1( entity, scalarLagrangeShapeFunctionSet );
  error = Dune::Fem::checkQuadratureConsistency( basisSet1, quadrature, true );
  if( error.two_norm() > eps )
  {
    std::cerr<<"Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy ): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );
  }

  error = 0;
  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLegendreShapeFunctionSetType >
  basisSet2( entity, scalarLegendreShapeFunctionSet );
  error = Dune::Fem::checkQuadratureConsistency( basisSet2, quadrature, true );
  if( error.two_norm() > eps )
  {
    std::cerr<<"Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy ): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LegendreShapeFunctionSet > test failed." );
  }

  error = 0;
  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarOrthonormalShapeFunctionSetType >
  basisSet3( entity, scalarOrthonormalShapeFunctionSet );
  error = Dune::Fem::checkQuadratureConsistency( basisSet3, quadrature, false );
  if( error.two_norm() > eps )
  {
    std::cerr<<"Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy ): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LegendreShapeFunctionSet > test failed." );
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

  if( gridPart.template begin< 0 >() == gridPart.template end< 0 >() )
    return 1;

  traverse< GridPartType, 1 >( gridPart );
  traverse< GridPartType, 2 >( gridPart );
  traverse< GridPartType, 3 >( gridPart );
  traverse< GridPartType, 4 >( gridPart );

  return 0;
}
