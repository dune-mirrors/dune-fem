#include <config.h>

#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/testgrid.hh>


// return random number in given interval
template< class FieldType >
FieldType random ( FieldType min, FieldType max )
{
  FieldType r = FieldType( rand() )/FieldType( RAND_MAX );
  return ( min + r*( max - min ) );
}


// compare basis function sets in all quadrature points
template< class BasisFunctionSetType, class VectorialBasisFunctionSetType, class QuadratureType >
void compareQuadrature ( const BasisFunctionSetType &basisFunctionSet, 
                         const VectorialBasisFunctionSetType &vectorialBasisFunctionSet,
                         const QuadratureType &quadrature )
{
  // check order of basis function set
  if( vectorialBasisFunctionSet.order() != basisFunctionSet.order() )
    DUNE_THROW( Dune::InvalidStateException, "order() does not match" );

  // check size of basis function set
  if( vectorialBasisFunctionSet.size() != basisFunctionSet.size() )
    DUNE_THROW( Dune::InvalidStateException, "size() does not match" );

  // get types
  typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  // init dof vectors
  const std::size_t size = basisFunctionSet.size();
  std::vector< RangeFieldType > dofs( basisFunctionSet.size() );
  for( std::size_t i = 0; i < size; ++i )
  {
    dofs[ i ] = random< RangeFieldType >( 1e-3, 1e3 );
    // std::cout << "dof[ " << i << " ] = " << dofs[ i ] << std::endl;
  }

  const std::size_t nop = quadrature.nop();
  for( std::size_t qp = 0; qp < nop; ++qp )
  {
    // check evaluate methods
    {
      RangeType value[ 2 ];
      basisFunctionSet.evaluateAll( quadrature[ qp ], dofs, value[ 0 ] );
      vectorialBasisFunctionSet.evaluateAll( quadrature[ qp ], dofs, value[ 1 ] );

      std::vector< RangeFieldType > d[ 2 ];
      d[ 0 ].resize( size, RangeFieldType( 0 ) );
      d[ 1 ].resize( size, RangeFieldType( 0 ) );
      basisFunctionSet.axpy( quadrature[ qp ], value[ 0 ], d[ 0 ] );
      vectorialBasisFunctionSet.axpy( quadrature[ qp ], value[ 0 ], d[ 1 ] );

      std::vector< RangeType > values[ 2 ];
      values[ 0 ].resize( size );
      values[ 1 ].resize( size );
      basisFunctionSet.evaluateAll( quadrature[ qp ], values[ 0 ] );
      vectorialBasisFunctionSet.evaluateAll( quadrature[ qp ], values[ 1 ] );
    }

    // check jacobian methods
    {
      JacobianRangeType jacobian[ 2 ];
      basisFunctionSet.jacobianAll( quadrature[ qp ], dofs, jacobian[ 0 ] );
      vectorialBasisFunctionSet.jacobianAll( quadrature[ qp ], dofs, jacobian[ 1 ] );

      std::vector< RangeFieldType > d[ 2 ];
      d[ 0 ].resize( size, RangeFieldType( 0 ) );
      d[ 1 ].resize( size, RangeFieldType( 0 ) );
      basisFunctionSet.axpy( quadrature[ qp ], jacobian[ 0 ], d[ 0 ] );
      vectorialBasisFunctionSet.axpy( quadrature[ qp ], jacobian[ 0 ], d[ 1 ] );

      std::vector< JacobianRangeType > jacobians[ 2 ];
      jacobians[ 0 ].resize( size );
      jacobians[ 1 ].resize( size );
      basisFunctionSet.jacobianAll( quadrature[ qp ], jacobians[ 0 ] );
      vectorialBasisFunctionSet.jacobianAll( quadrature[ qp ], jacobians[ 1 ] );
    }

    // check hessian methods
    {
      HessianRangeType hessian[ 2 ];
      basisFunctionSet.hessianAll( quadrature[ qp ], dofs, hessian[ 0 ] );
      vectorialBasisFunctionSet.hessianAll( quadrature[ qp ], dofs, hessian[ 1 ] );

      std::vector< HessianRangeType > hessians[ 2 ];
      hessians[ 0 ].resize( size );
      hessians[ 1 ].resize( size );
      basisFunctionSet.hessianAll( quadrature[ qp ], hessians[ 0 ] );
      vectorialBasisFunctionSet.hessianAll( quadrature[ qp ], hessians[ 1 ] );
    }
  }
}


template< class GridPartType >
void traverse ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimRange = DIMRANGE;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, dimRange > FunctionSpaceType;
  
  // create discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType space( gridPart );
  
  // create scalar discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< typename FunctionSpaceType::ScalarFunctionSpaceType, GridPartType, POLORDER > ScalarDiscreteFunctionSpaceType; 
  ScalarDiscreteFunctionSpaceType scalarSpace( gridPart );
  
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  const IteratorType end = space.end();
  for( IteratorType it = space.begin(); it != end; ++it )
  {
    // get entity
    const typename IteratorType::Entity &entity = *it;
    
    // get basis function set 
    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    const BasisFunctionSetType basisFunctionSet = space.basisFunctionSet( entity );
    
    // get scalar basis function set 
    typedef typename ScalarDiscreteFunctionSpaceType::BasisFunctionSetType ScalarBasisFunctionSetType;
    const ScalarBasisFunctionSetType scalarBasisFunctionSet = scalarSpace.basisFunctionSet( entity );
    
    // create vectorial basis function set
    typedef Dune::Fem::VectorialBasisFunctionSet< ScalarBasisFunctionSetType, typename BasisFunctionSetType::RangeType > VectorialBasisFunctionSetType;
    VectorialBasisFunctionSetType vectorialBasisFunctionSet( scalarBasisFunctionSet );
    
    // create quadrature
    assert( basisFunctionSet.order() == vectorialBasisFunctionSet.order() );
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
    QuadratureType quadrature( basisFunctionSet.entity(), basisFunctionSet.order() );

    // compare basis function sets in quadrature points
    compareQuadrature( basisFunctionSet, vectorialBasisFunctionSet, quadrature );
  }
}


int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  srand( time( NULL ) );

  typedef Dune::GridSelector::GridType GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();
  
  int refCount = Dune::Fem::Parameter::getValue< int >( "startLevel", 0 );
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( grid, refCount*refineStepsForHalf );
  
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );
  traverse( gridPart );
  
  return 0;
}
