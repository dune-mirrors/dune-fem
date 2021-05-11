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

  // type of error
  typedef Dune::FieldVector< double, 7 > ErrorType;

  // prepare shapefunctions
  ScalarLagrangeShapeFunctionSetType scalarLagrangeShapeFunctionSet( entity.type() );
  ScalarLegendreShapeFunctionSetType scalarLegendreShapeFunctionSet( polorder );

  // CXXFLAGS="-std=c++17 -O3 -ffast-math -march=native" leads to a small
  // increase in the error so that (with g++-9.3) fails with eps=1e-6
  double eps = 2e-6;

  ErrorType error( 0 );
  typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLagrangeShapeFunctionSetType > ScalarBasisFunctionSetType1;
  typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLegendreShapeFunctionSetType > ScalarBasisFunctionSetType2;
  typedef Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType2, typename FunctionSpaceType::RangeType,
                                                Dune::Fem::HorizontalDofAlignment > BasisFunctionSetType3;

  ScalarBasisFunctionSetType1 scalarBasisSet1( entity, scalarLagrangeShapeFunctionSet );
  ScalarBasisFunctionSetType2 scalarBasisSet2( entity, scalarLegendreShapeFunctionSet );
  BasisFunctionSetType3 basisFunctionSet3( scalarBasisSet2 );

  Dune::Fem::TupleBasisFunctionSet< Dune::Fem::TupleSpaceProduct,
              ScalarBasisFunctionSetType1, ScalarBasisFunctionSetType2, BasisFunctionSetType3 >
  basisSet( scalarBasisSet1, scalarBasisSet2, basisFunctionSet3 );

  error = Dune::Fem::checkQuadratureConsistency( basisSet, quadrature, true );
  if( error.infinity_norm() > eps )
  {
    std::cerr<<"Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " TupleBasisFunctionSet< ShapeFunctionSet1, ShapeFunctionSet2, ShapeFunctionSet3 > test failed." );
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
