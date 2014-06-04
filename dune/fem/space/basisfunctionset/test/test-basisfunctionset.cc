#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/testgrid.hh>


// random
// ------

template< class FieldType >
FieldType random ( FieldType a, FieldType b )
{
  FieldType delta = FieldType( rand() )/FieldType( RAND_MAX );
  return ( a + delta*( b - a ) );
}


// compareQuadrature
// -----------------

template< class BasisFunctionSetType, class QuadratureType >
Dune::FieldVector< typename BasisFunctionSetType::RangeType::value_type, 3 >
  checkQuadratureConsistency ( const BasisFunctionSetType &basisFunctionSet, 
                               const QuadratureType &quadrature )
{
  // get types
  typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const int dimRange = FunctionSpaceType :: dimRange;

  const std::size_t size = basisFunctionSet.size();

  // init random dof vectors
  std::vector< RangeFieldType > dofs( size );
  for( std::size_t i = 0; i < size; ++i )
    dofs[ i ] = random< RangeFieldType >( 1e-3, 1e3 );

  // return value
  Dune::FieldVector< RangeFieldType, 3 > ret;

  const std::size_t nop = quadrature.nop();
  for( std::size_t qp = 0; qp < nop; ++qp )
  {
    // check evaluate methods
    {
      RangeFieldType error = 0;
      RangeType a,b;
      std::vector< RangeType > values;
      values.resize( basisFunctionSet.size() );

      basisFunctionSet.evaluateAll( quadrature[ qp ], dofs, a );
      basisFunctionSet.evaluateAll( quadrature[ qp ], values );

      for( int i =0; i< values.size(); ++i )
        b.axpy( dofs[ i ], values[ i ] );
      a -= b;

      error = a.two_norm();
      ret[ 0 ] = std::max( ret[ 0 ], std::abs( error ) );
    }

    // check jacobian methods
    {
      RangeFieldType error = 0;
      JacobianRangeType a,b;
      std::vector< JacobianRangeType > values;
      values.resize( basisFunctionSet.size() );

      basisFunctionSet.jacobianAll( quadrature[ qp ], dofs, a );
      basisFunctionSet.jacobianAll( quadrature[ qp ], values );

      for( int i =0; i< values.size(); ++i )
        b.axpy( dofs[ i ], values[ i ] );
      a -= b;

      error = a.frobenius_norm();
      ret[ 1 ] = std::max( ret[ 1 ], error );
    }

    // check hessian methods
    {
      RangeFieldType error = 0;
      HessianRangeType a,b;
      std::vector< HessianRangeType > values;
      values.resize( basisFunctionSet.size() );

      basisFunctionSet.hessianAll( quadrature[ qp ], dofs, a );
      basisFunctionSet.hessianAll( quadrature[ qp ], values );


      for( int r =0; r < dimRange; ++r )
      {
        for( int i =0; i< values.size(); ++i )
          b[ r ].axpy( dofs[ i ], values[ i ][ r ] );
      }

      for( int r =0; r < dimRange; ++r )
        a[ r ] -= b[ r ];

      for( int r =0; r < dimRange; ++r )
        error += a[ r ].frobenius_norm2();

      ret[ 2 ] = std::max( ret[ 1 ], std::sqrt( error ) );
    }
  }

  return ret;
}


template< class GridPartType >
void traverse ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimRange = DIMRANGE;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, dimRange > FunctionSpaceType;
  
  // create discrete function space, needs to be changed
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType space( gridPart );
  
  typedef Dune::FieldVector< typename FunctionSpaceType::RangeFieldType, 3 > ErrorType;
  ErrorType error( 0 );
  
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  const IteratorType end = space.end();
  for( IteratorType it = space.begin(); it != end; ++it )
  {
    // get entity
    const typename IteratorType::Entity &entity = *it;
    
    // get basis function set 
    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    typedef typename DiscreteFunctionSpaceType::ShapeFunctionSetType ShapeFunctionSetType;
    const BasisFunctionSetType basisFunctionSet = space.basisFunctionSet( entity );
    
    // create quadrature
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
    QuadratureType quadrature( basisFunctionSet.entity(), basisFunctionSet.order() );

    ErrorType localError;
    {
      // compare basis function sets in quadrature points
      Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > defaultBasisFunctionSet( entity,
          basisFunctionSet.shapeFunctionSet() );
      localError = checkQuadratureConsistency( defaultBasisFunctionSet, quadrature );

      for( int i = 0; i < 3; ++i )
        error[ i ] = std::max( error[ i ], localError[ i ] );
    }

    // update error
    for( int i = 0; i < 3; ++i )
      error[ i ] = std::max( error[ i ], localError[ i ] );
  }

  // print error
  std::cout << "Error (evaluateAll, jacobianAll, hessianAll ) = " << error << std::endl;
}


int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( argc >= 2 ? argv[ 1 ] : "parameter" );

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
