#ifndef DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH
#define DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH

#include <algorithm>
#include <functional>
#include <iostream>
#include <random>

#include <dune/common/fvector.hh>
#include <dune/fem/common/coordinate.hh>

namespace Dune
{

  namespace Fem
  {

    // checkQuadratureConsistency
    // --------------------------

    template< class BasisFunctionSet, class Quadrature >
    Dune::FieldVector< typename BasisFunctionSet::RangeType::value_type, 7 >
    checkQuadratureConsistency ( const BasisFunctionSet &basisFunctionSet, const Quadrature &quadrature,
                                 bool supportsHessians = false )
    {
      // get types
      typedef typename BasisFunctionSet::FunctionSpaceType FunctionSpaceType;
      typedef typename BasisFunctionSet::DomainType DomainType;
      typedef typename BasisFunctionSet::RangeType RangeType;
      typedef typename BasisFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSet::HessianRangeType HessianRangeType;
      typedef typename BasisFunctionSet::EntityType EntityType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      static const int dimRange = FunctionSpaceType::dimRange;

      const EntityType &entity = basisFunctionSet.entity();
      const auto refElement = basisFunctionSet.referenceElement();

      if( entity.type() != refElement.type() )
        DUNE_THROW( Dune::InvalidStateException, "GeometryType of referenceElement and entity mismatch for this basisFunctionSet" );

      int order = basisFunctionSet.order();
      // prevent warning about unused params
      (void) order;

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

      HessianRangeType hessianFactor;
      for( std::size_t j = 0; j < JacobianRangeType::rows; ++j )
        for( std::size_t k = 0; k < JacobianRangeType::cols; ++k )
          for( std::size_t l = 0; l < JacobianRangeType::cols; ++l )
            hessianFactor[ j ][ k ][ l ] = random();

      // return value
      Dune::FieldVector< RangeFieldType, 7 > ret( 0 );

      const std::size_t nop = quadrature.nop();

      std::vector< RangeType > aVec( nop );
      std::vector< JacobianRangeType > aJac( nop );

      basisFunctionSet.evaluateAll( quadrature, dofs, aVec );
      basisFunctionSet.jacobianAll( quadrature, dofs, aJac );

      for( std::size_t qp = 0; qp < nop; ++qp )
      {
        // check evaluate methods
        {
          RangeType  b;
          RangeType& a = aVec[ qp ];

          std::vector< RangeType > values( basisFunctionSet.size() );
          basisFunctionSet.evaluateAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
          {
            b.axpy( dofs[ i ], values[ i ] );
          }

          ret[ 0 ] = std::max( ret[ 0 ], (a - b).two_norm() );
        }

        // check jacobian methods
        {
          JacobianRangeType b;
          JacobianRangeType& a = aJac[ qp ];

          std::vector< JacobianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.jacobianAll( quadrature[ qp ], values );
          for( std::size_t i = 0; i < values.size(); ++i )
            b.axpy( dofs[ i ], values[ i ] );

          a -= b;
          ret[ 1 ] = std::max( ret[ 1 ], a.frobenius_norm() );
        }

        // check hessian methods
        if( supportsHessians )
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
          DomainType x = coordinate( quadrature[ qp ] );
          std::vector< RangeFieldType > r1( dofs );
          std::vector< RangeFieldType > r2( dofs );

          basisFunctionSet.axpy( x, jacobianFactor, r1 );

          std::vector< JacobianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.jacobianAll( x, values );
          for( std::size_t i = 0; i < values.size(); ++i )
            for( int j = 0; j < JacobianRangeType::rows; ++j )
              r2[ i ] += jacobianFactor[ j ] * values[ i ][ j ];

          RangeFieldType error = 0;
          for( std::size_t i = 0; i < values.size(); ++i )
            error += std::abs( r2[ i ] - r1[ i ] );

          ret[ 4 ] = std::max( ret[ 4 ], error );
        }


        // check hessian axpy method
        if( supportsHessians )
        {
          DomainType x = coordinate( quadrature[ qp ] );
          std::vector< RangeFieldType > r1( dofs );
          std::vector< RangeFieldType > r2( dofs );

          basisFunctionSet.axpy( x, hessianFactor, r1 );

          std::vector< HessianRangeType > values( basisFunctionSet.size() );
          basisFunctionSet.hessianAll( x, values );
          for( std::size_t i = 0; i < values.size(); ++i )
            for( int j = 0; j < JacobianRangeType::rows; ++j )
              for( std::size_t k = 0; k < JacobianRangeType::cols; ++k )
                for( std::size_t l = 0; l < JacobianRangeType::cols; ++l )
                  r2[ i ] += hessianFactor[ j ][ k ][ l ] * values[ i ][ j ][ k ][ l ];

          RangeFieldType error = 0;
          for( std::size_t i = 0; i < values.size(); ++i )
            error += std::abs( r2[ i ] - r1[ i ] );

          ret[ 5 ] = std::max( ret[ 5 ], error );
        }

        // check value, jaacobian axpy method
        {
          DomainType x = coordinate( quadrature[ qp ] );
          std::vector< RangeFieldType > r1( dofs );
          std::vector< RangeFieldType > r2( dofs );

          basisFunctionSet.axpy( x, valueFactor, jacobianFactor, r1 );

          std::vector< RangeType > values( basisFunctionSet.size() );
          std::vector< JacobianRangeType > jacobians( basisFunctionSet.size() );

          basisFunctionSet.evaluateAll( x, values );
          basisFunctionSet.jacobianAll( x, jacobians );

          for( std::size_t i = 0; i < values.size(); ++i )
          {
            r2[ i ] += valueFactor * values[ i ];
            for( int j = 0; j < JacobianRangeType::rows; ++j )
              r2[ i ] += jacobianFactor[ j ] * jacobians[ i ][ j ];
          }

          RangeFieldType error = 0;
          for( std::size_t i = 0; i < values.size(); ++i )
            error += std::abs( r2[ i ] - r1[ i ] );

          ret[ 6 ] = std::max( ret[ 6 ], error );
        }
      }

      return ret;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_TEST_CHECKBASISFUNCTIONSET_HH
