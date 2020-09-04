#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#ifdef USE_BASEFUNCTIONSET_CODEGEN
#include <dune/fem/space/basisfunctionset/default_codegen.hh>
#endif

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/interpolationquadrature.hh>

#include <dune/fem/space/shapefunctionset/lagrange.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>
#include <dune/fem/space/shapefunctionset/caching.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/simple.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include <dune/fem/space/lagrange.hh>

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
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0, Dune::Fem::DefaultQuadratureTraits > QuadratureType;
  QuadratureType quadrature( entity, polorder );
#if HAVE_DUNE_LOCALFUNCTIONS
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0, Dune::Fem::GaussLobattoQuadratureTraits > GLQuadratureType;
  GLQuadratureType glQuad( entity, polorder );
  GLQuadratureType gl2Quad( entity, polorder+1 );
#endif

  // needs a geometry type to construct
  typedef Dune::Fem::CachingShapeFunctionSet < Dune::Fem::LagrangeShapeFunctionSet< ScalarFunctionSpaceType, polorder > > ScalarLagrangeShapeFunctionSetType;

  // needs a geometry type to construct
  typedef Dune::Fem::CachingShapeFunctionSet < Dune::Fem::LagrangeShapeFunctionSet< ScalarFunctionSpaceType, polorder > > ScalarLagrangeShapeFunctionSetType;

  // needs an order to construct
  typedef Dune::Fem::CachingShapeFunctionSet< Dune::Fem::LegendreShapeFunctionSet< ScalarFunctionSpaceType > > ScalarLegendreShapeFunctionSetType;

  // needs a geometry type to construct
  typedef Dune::Fem::CachingShapeFunctionSet< Dune::Fem::OrthonormalShapeFunctionSet< ScalarFunctionSpaceType > > ScalarOrthonormalShapeFunctionSetType;

  // type of error
  typedef Dune::FieldVector< double, 7 > ErrorType;

  // prepare shapefunctions
  typename ScalarLagrangeShapeFunctionSetType::ShapeFunctionSetType lagset( entity.type() );
  ScalarLagrangeShapeFunctionSetType scalarLagrangeShapeFunctionSet( entity.type(), lagset );

  typename ScalarLegendreShapeFunctionSetType::ShapeFunctionSetType implset( polorder);
  ScalarLegendreShapeFunctionSetType scalarLegendreShapeFunctionSet(
      entity.type(), implset );

  typename ScalarOrthonormalShapeFunctionSetType::ShapeFunctionSetType orthoimplset(entity.type(), polorder);
  ScalarOrthonormalShapeFunctionSetType scalarOrthonormalShapeFunctionSet( entity.type(), orthoimplset );

  double eps = 1e-7;

  ErrorType error( 0 );

#if HAVE_DUNE_LOCALFUNCTIONS && not defined(USE_BASEFUNCTIONSET_CODEGEN)
  typedef Dune::Fem::LagrangeFiniteElementMap< ScalarFunctionSpaceType, GridPartType, Dune::GaussLobattoPointSet > LFEMap;
  typedef Dune::Fem::LocalFunctionsShapeFunctionSet< typename LFEMap::LocalFiniteElementType::Traits::LocalBasisType, LFEMap::pointSetId > LocalFunctionsShapeFunctionSetType;
  typedef Dune::Fem::CachingShapeFunctionSet < LocalFunctionsShapeFunctionSetType > ScalarLocalFiniteElementShapeFunctionSetType;

  // prepare shapefunctions
  LFEMap lfemap( gridPart, polorder );

  LocalFunctionsShapeFunctionSetType lfeSet( std::get< 1 >(lfemap( entity )) );
  ScalarLocalFiniteElementShapeFunctionSetType scalarLocalFiniteElementFunctionSet( entity.type(), lfeSet );

  // default basis function set
  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLocalFiniteElementShapeFunctionSetType >
  basisSet0( entity, scalarLocalFiniteElementFunctionSet );
  //std::cout << "Check basisSet0 " << std::endl;

  {
    error = Dune::Fem::checkQuadratureConsistency( basisSet0, glQuad, false );
    if( error.infinity_norm() > eps )
    {
      std::cerr<<"set1: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
      DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );
    }
  }

  {
    error = Dune::Fem::checkQuadratureConsistency( basisSet0, gl2Quad, false );
    if( error.infinity_norm() > eps )
    {
      std::cerr<<"set1: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
      DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );
    }
  }

  error = Dune::Fem::checkQuadratureConsistency( basisSet0, quadrature, false );
  if( error.infinity_norm() > eps )
  {
    std::cerr<<"set1: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );
  }

#endif

  // default basis function set
  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLagrangeShapeFunctionSetType >
  basisSet1( entity, scalarLagrangeShapeFunctionSet );
  //std::cout << "Check basisSet1 " << std::endl;
  error = Dune::Fem::checkQuadratureConsistency( basisSet1, quadrature, false );
  if( error.infinity_norm() > eps )
  {
    std::cerr<<"set1: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LagrangeShapeFunctionSet > test failed." );
  }

  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarLegendreShapeFunctionSetType >
  basisSet2( entity, scalarLegendreShapeFunctionSet );
  //std::cout << "Check basisSet2 " << std::endl;
  error = Dune::Fem::checkQuadratureConsistency( basisSet2, quadrature, false );
  if( error.infinity_norm() > eps )
  {
    std::cerr<<"set2: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
    DUNE_THROW( Dune::InvalidStateException, " DefaultBasisFunctionSet< LegendreShapeFunctionSet > test failed." );
  }

  Dune::Fem::DefaultBasisFunctionSet< EntityType, ScalarOrthonormalShapeFunctionSetType >
  basisSet3( entity, scalarOrthonormalShapeFunctionSet );
  //std::cout << "Check basisSet3 " << std::endl;
  error = Dune::Fem::checkQuadratureConsistency( basisSet3, quadrature, false );
  if( error.infinity_norm() > eps )
  {
    std::cerr<<"set3: Errors( evaluate, jacobian, hessian, value axpy, jacobian axpy, hessian axpy, v+j axpy): "<< error <<std::endl;
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

  if( gridPart.begin< 0 >() == gridPart.end< 0 >() )
    return 1;

  traverse< GridPartType, 1 >( gridPart );
  traverse< GridPartType, 2 >( gridPart );
  traverse< GridPartType, 3 >( gridPart );
  traverse< GridPartType, 4 >( gridPart );

  return 0;
}
