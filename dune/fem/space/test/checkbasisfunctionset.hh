#ifndef DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH
#define DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH

#include <algorithm>
#include <iostream>
#include <random>

#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {

    // checkQuadratureConsistency
    // --------------------------

    template< class BasisFunctionSet, class Quadrature >
    Dune::FieldVector< typename BasisFunctionSet::RangeType::value_type, 5 >
    checkQuadratureConsistency ( const BasisFunctionSet &basisFunctionSet, const Quadrature &quadrature,
                                 bool supporsHessians = false )
    {
      // get types
      typedef typename BasisFunctionSet::FunctionSpaceType FunctionSpaceType;
      typedef typename BasisFunctionSet::DomainType DomainType;
      typedef typename BasisFunctionSet::RangeType RangeType;
      typedef typename BasisFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSet::HessianRangeType HessianRangeType;
      typedef typename BasisFunctionSet::ReferenceElementType ReferenceElementType;
      typedef typename BasisFunctionSet::EntityType EntityType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      static const int dimRange = FunctionSpaceType::dimRange;

      const ReferenceElementType &refElement = basisFunctionSet.referenceElement();

      int order = basisFunctionSet.order();

      const std::size_t size = basisFunctionSet.size();

      // init random dof vectors
      std::uniform_real_distribution< RangeFieldType > distribution( 1e-3, 1e3 );
      std::default_random_engine randomEngine;
      auto random = std::bind( distribution, randomEngine );

      std::vector< RangeFieldType > dofs( size );
      for( RangeFieldType &dof : dofs )
        dof = random();

      RangeType valueFactor;
      for( RangeFieldType &v : valueFactor )
        v = random();

      JacobianRangeType jacobianFactor;
      for( std::size_t j = 0; j < JacobianRangeType::rows; ++j )
        for( std::size_t k = 0; k < JacobianRangeType::cols; ++k )
          jacobianFactor[ j ][ k ] = random();

      // return value
      Dune::FieldVector< RangeFieldType, 5 > ret;

      const std::size_t nop = quadrature.nop();
      for( std::size_t qp = 0; qp < nop; ++qp )
      {
        // check evaluate methods
        {
          RangeType a, b;

          basisFunctionSet.evaluateAll( quadrature[ qp ], dofs, a );

          std::vector< RangeType > values( basisFunctionSet.size() );
          basisFunctionSet.evaluateAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
            b.axpy( dofs[ i ], values[ i ] );

          ret[ 0 ] = std::max( ret[ 0 ], (a - b).two_norm() );
        }

        // check jacobian methods
        {
          JacobianRangeType a, b;

          basisFunctionSet.jacobianAll( quadrature[ qp ], dofs, a );

          std::vector< JacobianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.jacobianAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
            b.axpy( dofs[ i ], values[ i ] );

          a -= b;
          ret[ 1 ] = std::max( ret[ 1 ], a.frobenius_norm() );
        }

        // check hessian methods
        if( supporsHessians )
        {
          HessianRangeType a, b;

          basisFunctionSet.hessianAll( quadrature[ qp ], dofs, a );

          std::vector< HessianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.hessianAll( quadrature[ qp ], values );
          for( int r = 0; r < dimRange; ++r )
            for( std::size_t i = 0; i < values.size(); ++i )
              b[ r ].axpy( dofs[ i ], values[ i ][ r ] );

          RangeFieldType error( 0 );
          for( int r = 0; r < dimRange; ++r )
            a[ r ] -= b[ r ];
          for( int r = 0; r < dimRange; ++r )
            error += a[ r ].frobenius_norm2();

          ret[ 2 ] = std::max( ret[ 2 ], std::sqrt( error ) );
        }

        // check value axpy method
        {
          std::vector< RangeFieldType > r1( dofs );
          std::vector< RangeFieldType > r2( dofs );

          basisFunctionSet.axpy( quadrature[ qp ], valueFactor, r1 );

          std::vector< RangeType > values( basisFunctionSet.size() );
          basisFunctionSet.evaluateAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
            r2[ i ] += valueFactor * values[ i ];

          RangeFieldType error = 0;
          for( std::size_t i = 0; i < values.size(); ++i )
            error += std::abs( r2[ i ] - r1[ i ] );

          ret[ 3 ] = std::max( ret[ 3 ], error );
        }

        // check jacobian axpy method
        {
          std::vector< RangeFieldType > r1( dofs );
          std::vector< RangeFieldType > r2( dofs );

          basisFunctionSet.axpy( quadrature[ qp ], jacobianFactor, r1 );

          std::vector< JacobianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.jacobianAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
            for( int j = 0; j < JacobianRangeType::rows; ++j )
              r2[ i ] += jacobianFactor[ j ] * values[ i ][ j ];

          RangeFieldType error = 0;
          for( std::size_t i = 0; i < values.size(); ++i )
            error += std::abs( r2[ i ] - r1[ i ] );

          ret[ 4 ] = std::max( ret[ 4 ], error );
        }
      }

      return ret;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH
