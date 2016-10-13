#ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
#define DUNE_FEM_SCHEMES_GALERKIN_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/ftraits.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>

#include <dune/fem/schemes/solver.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      // GalerkinOperator
      // ----------------

      template< class DiscreteFunctionSpace >
      struct GalerkinOperator
      {
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
        typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

        typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

      private:
        typedef CachingQuadrature< GridPartType, 0 > InteriorQuadratureType;
        typedef QuadraturePointWrapper< InteriorQuadratureType > InteriorQuadraturePointType;

        typedef CachingQuadrature< GridPartType, 1 > SurfaceQuadratureType;
        typedef QuadraturePointWrapper< SurfaceQuadratureType > SurfaceQuadraturePointType;

        typedef std::tuple< RangeType, JacobianRangeType > ValueType;
        typedef std::tuple< std::vector< RangeType >, std::vector< JacobianRangeType > > ValueVectorType;

        template< class LocalFunction, class Point >
        static ValueType value ( const LocalFunction &u, const Point &x )
        {
          ValueType value;
          u.evaluate( x, std::get< 0 >( value ) );
          u.jacobian( x, std::get< 1 >( value ) );
          return value;
        }

        template< class Basis, class Point >
        static void values ( const Basis &basis, const Point &x, ValueVectorType &phi )
        {
          basis.evaluateAll( x, std::get< 0 >( phi ) );
          basis.jacobianAll( x, std::get< 1 >( phi ) );
        }

        // interior integrand

      private:
        template< class Integrands >
        static std::true_type chkHasInteriorIntegrand ( const Integrands &, decltype( std::declval< const Integrands & >().interior( std::declval< const InteriorQuadraturePointType & >(), std::declval< const ValueType & >() ) ) * = nullptr );

        static std::false_type chkHasInteriorIntegrand ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasInterior() ), bool >::value, int > = 0 >
        static std::true_type chkHasInteriorIntegrandCheck ( const Integrands & );

        static std::false_type chkHasInteriorIntegrandCheck ( ... );

      public:
        template< class Integrands >
        struct HasInteriorIntegrand
          : public decltype( chkHasInteriorIntegrand( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands >
        struct HasInteriorIntegrandCheck
          : public decltype( chkHasInteriorIntegrandCheck( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands, std::enable_if_t< HasInteriorIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasInteriorIntegrand ( const Integrands &integrands )
        {
          return integrands.hasInterior();
        }

        template< class Integrands, std::enable_if_t< !HasInteriorIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasInteriorIntegrand ( const Integrands &integrands )
        {
          return HasInteriorIntegrand< Integrands >::value;
        }

        template< class Integrands, class W, std::enable_if_t< HasInteriorIntegrand< Integrands >::value, int > = 0 >
        static void addInteriorIntegrand ( const Integrands &integrands, const InteriorQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, W &w )
        {
          ValueType integrand = integrands.interior( qp, u );
          std::get< 0 >( integrand ) *= weight;
          std::get< 1 >( integrand ) *= weight;
          w.axpy( qp, std::get< 0 >( integrand ), std::get< 1 >( integrand ) );
        }

        template< class Integrands, class W, std::enable_if_t< !HasInteriorIntegrand< Integrands >::value, int > = 0 >
        static void addInteriorIntegrand ( const Integrands &integrands, const InteriorQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, W &w )
        {}

        template< class Integrands, class J, std::enable_if_t< HasInteriorIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedInteriorIntegrand ( const Integrands &integrands, const InteriorQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, const ValueVectorType &phi, std::size_t cols, J &j )
        {
          auto integrand = integrands.linearizedInterior( qp, u );
          for( std::size_t col = 0; col < cols; ++col )
          {
            ValueType intPhi = integrand( std::make_tuple( std::get< 0 >( phi )[ col ], std::get< 1 >( phi )[ col ] ) );
            j.column( col ).axpy( std::get< 0 >( phi ), std::get< 1 >( phi ), std::get< 0 >( intPhi ), std::get< 1 >( intPhi ), weight );
          }
        }

        template< class Integrands, class J, std::enable_if_t< !HasInteriorIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedInteriorIntegrand ( const Integrands &integrands, const InteriorQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, const ValueVectorType &phi, std::size_t cols, J &j )
        {}

        // boundary integrand

      private:
        template< class Integrands >
        static std::true_type chkHasBoundaryIntegrand ( const Integrands &, decltype( std::declval< const Integrands & >().boundary( std::declval< const SurfaceQuadraturePointType & >(), std::declval< const ValueType & >() ) ) * = nullptr );

        static std::false_type chkHasBoundaryIntegrand ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasBoundary() ), bool >::value, int > = 0 >
        static std::true_type chkHasBoundaryIntegrandCheck ( const Integrands & );

        static std::false_type chkHasBoundaryIntegrandCheck ( ... );

      public:
        template< class Integrands >
        struct HasBoundaryIntegrand
          : public decltype( chkHasBoundaryIntegrand( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands >
        struct HasBoundaryIntegrandCheck
          : public decltype( chkHasBoundaryIntegrandCheck( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands, std::enable_if_t< HasBoundaryIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasBoundaryIntegrand ( const Integrands &integrands )
        {
          return integrands.hasBoundary();
        }

        template< class Integrands, std::enable_if_t< !HasBoundaryIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasBoundaryIntegrand ( const Integrands &integrands )
        {
          return HasBoundaryIntegrand< Integrands >::value;
        }

        template< class Integrands, class W, std::enable_if_t< HasBoundaryIntegrand< Integrands >::value, int > = 0 >
        static void addBoundaryIntegrand ( const Integrands &integrands, const SurfaceQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, W &w )
        {
          ValueType integrand = integrands.boundary( qp, u );
          std::get< 0 >( integrand ) *= weight;
          std::get< 1 >( integrand ) *= weight;
          w.axpy( qp, std::get< 0 >( integrand ), std::get< 1 >( integrand ) );
        }

        template< class Integrands, class W, std::enable_if_t< !HasBoundaryIntegrand< Integrands >::value, int > = 0 >
        static void addBoundaryIntegrand ( const Integrands &integrands, const SurfaceQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, W &w )
        {}

        template< class Integrands, class J, std::enable_if_t< HasBoundaryIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedBoundaryIntegrand ( const Integrands &integrands, const SurfaceQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, const ValueVectorType &phi, std::size_t cols, J &j )
        {
          auto integrand = integrands.linearizedBoundary( qp, u );
          for( std::size_t col = 0; col < cols; ++col )
          {
            ValueType intPhi = integrand( std::make_tuple( std::get< 0 >( phi )[ col ], std::get< 1 >( phi )[ col ] ) );
            j.column( col ).axpy( std::get< 0 >( phi ), std::get< 1 >( phi ), std::get< 0 >( intPhi ), std::get< 1 >( intPhi ), weight );
          }
        }

        template< class Integrands, class J, std::enable_if_t< !HasBoundaryIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedBoundaryIntegrand ( const Integrands &integrands, const SurfaceQuadraturePointType &qp, const RangeFieldType &weight, const ValueType &u, const ValueVectorType &phi, std::size_t cols, J &j )
        {}

        // skeleton integrand

      private:
        template< class Integrands >
        static std::true_type chkHasSkeletonIntegrand ( const Integrands &, decltype( std::declval< const Integrands & >().skeleton( std::declval< const SurfaceQuadraturePointType & >(), std::declval< const ValueType & >(), std::declval< const SurfaceQuadraturePointType & >(), std::declval< const ValueType & >() ) ) * = nullptr );

        static std::false_type chkHasSkeletonIntegrand ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasSkeleton() ), bool >::value, int > = 0 >
        static std::true_type chkHasSkeletonIntegrandCheck ( const Integrands & );

        static std::false_type chkHasSkeletonIntegrandCheck ( ... );

      public:
        template< class Integrands >
        struct HasSkeletonIntegrand
          : public decltype( chkHasSkeletonIntegrand( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands >
        struct HasSkeletonIntegrandCheck
          : public decltype( chkHasSkeletonIntegrandCheck( std::declval< const Integrands & >() ) )
        {};

        template< class Integrands, std::enable_if_t< HasSkeletonIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasSkeletonIntegrand ( const Integrands &integrands )
        {
          return integrands.hasSkeleton();
        }

        template< class Integrands, std::enable_if_t< !HasSkeletonIntegrandCheck< Integrands >::value, int > = 0 >
        static bool hasSkeletonIntegrand ( const Integrands &integrands )
        {
          return HasSkeletonIntegrand< Integrands >::value;
        }

        template< class Integrands, class QP, class W, std::enable_if_t< HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, W &wIn )
        {
          std::pair< ValueType, ValueType > integrand = integrands.skeleton( qpIn, uIn, uOut );
          std::get< 0 >( integrand.first ) *= weight;
          std::get< 1 >( integrand.first ) *= weight;
          wIn.axpy( qpIn, std::get< 0 >( integrand.first ), std::get< 1 >( integrand.first ) );
        }

        template< class Integrands, class QP, class W, std::enable_if_t< HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, W &wIn, W &wOut )
        {
          std::pair< ValueType, ValueType > integrand = integrands.skeleton( qpIn, uIn, qpOut, uOut );
          std::get< 0 >( integrand.first ) *= weight;
          std::get< 1 >( integrand.first ) *= weight;
          wIn.axpy( qpIn, std::get< 0 >( integrand.first ), std::get< 1 >( integrand.first ) );
          std::get< 0 >( integrand.second ) *= weight;
          std::get< 1 >( integrand.second ) *= weight;
          wOut.axpy( qpOut, std::get< 0 >( integrand.second ), std::get< 1 >( integrand.second ) );
        }

        template< class Integrands, class QP, class W, std::enable_if_t< !HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, W &wIn )
        {}

        template< class Integrands, class QP, class W, std::enable_if_t< !HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, W &wIn, W &wOut )
        {}

        template< class Integrands, class QP, class J, std::enable_if_t< HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn )
        {
          auto integrand = integrands.skeleton( qpIn, uIn, qpOut, uOut );
          for( std::size_t col = 0; col < colsIn; ++col )
          {
            std::pair< ValueType, ValueType > intPhi = integrand.first( std::make_tuple( std::get< 0 >( phiIn )[ col ], std::get< 1 >( phiIn )[ col ] ) );
            jInIn.column( col ).axpy( std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), std::get< 0 >( intPhi.first ), std::get< 1 >( intPhi.first ), weight );
          }
          for( std::size_t col = 0; col < colsOut; ++col )
          {
            std::pair< ValueType, ValueType > intPhi = integrand.second( std::make_tuple( std::get< 0 >( phiOut )[ col ], std::get< 1 >( phiOut )[ col ] ) );
            jOutIn.column( col ).axpy( std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), std::get< 0 >( intPhi.first ), std::get< 1 >( intPhi.first ), weight );
          }
        }

        template< class Integrands, class QP, class J, std::enable_if_t< HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn, J &jInOut, J &jOutOut )
        {
          auto integrand = integrands.linearizedSkeleton( qpIn, uIn, uOut );
          for( std::size_t col = 0; col < colsIn; ++col )
          {
            std::pair< ValueType, ValueType > intPhi = integrand.first( std::make_tuple( std::get< 0 >( phiIn )[ col ], std::get< 1 >( phiIn )[ col ] ) );
            jInIn.column( col ).axpy( std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), std::get< 0 >( intPhi.first ), std::get< 1 >( intPhi.first ), weight );
            jInOut.column( col ).axpy( std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), std::get< 0 >( intPhi.second ), std::get< 1 >( intPhi.second ), weight );
          }
          for( std::size_t col = 0; col < colsOut; ++col )
          {
            std::pair< ValueType, ValueType > intPhi = integrand.second( std::make_tuple( std::get< 0 >( phiOut )[ col ], std::get< 1 >( phiOut )[ col ] ) );
            jOutIn.column( col ).axpy( std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), std::get< 0 >( intPhi.first ), std::get< 1 >( intPhi.first ), weight );
            jOutOut.column( col ).axpy( std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), std::get< 0 >( intPhi.second ), std::get< 1 >( intPhi.second ), weight );
          }
        }

        template< class Integrands, class QP, class J, std::enable_if_t< !HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn )
        {}

        template< class Integrands, class QP, class J, std::enable_if_t< !HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        static void addLinearizedSkeletonIntegrand ( const Integrands &integrands, const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn, J &jInOut, J &jOutOut )
        {}

        // interior integral

        template< class Integrands, class U, class W >
        void addInteriorIntegral ( Integrands &integrands, const U &u, W &w ) const
        {
          if( !integrands.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();
          for( const InteriorQuadraturePointType qp : InteriorQuadratureType( u.entity(), 2*w.order() ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.position() );
            addInteriorIntegrand( integrands, qp, weight, value( u, qp ), w );
          }
        }

        template< class Integrands, class U, class J >
        void addLinearizedInteriorIntegral ( Integrands &integrands, const U &u, ValueVectorType &phi, J &j ) const
        {
          if( !integrands.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();
          const auto &basis = discreteFunctionSpace().basisFunctionSet( u.entity() );
          for( const auto qp : InteriorQuadratureType( u.entity(), 2*basis.order() ) )
          {
            const auto weight = qp.weight() * geometry.integrationElement( qp.position() );
            values( basis, qp, phi );
            addLinearizedInteriorIntegrand( integrands, qp, weight, value( u, qp ), phi, basis.size(), j );
          }
        }

        // boundary integral

        template< class Integrands, class Intersection, class U, class W >
        void addBoundaryIntegral ( Integrands &integrands, const Intersection &intersection, const U &u, W &w ) const
        {
          if( !integrands.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          for( const SurfaceQuadraturePointType qp : SurfaceQuadratureType( gridPart(), intersection, 2*w.order(), SurfaceQuadratureType::INSIDE ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.localPosition() );
            addBoundaryIntegrand( integrands, qp, weight, value( u, qp ), w );
          }
        }

        template< class Integrands, class Intersection, class U, class J >
        void addLinearizedBoundaryIntegral ( Integrands &integrands, const Intersection &intersection, const U &u, ValueVectorType &phi, J &j ) const
        {
          if( !integrands.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          const auto &basis = discreteFunctionSpace().basisFunctionSet( u.entity() );
          for( const SurfaceQuadraturePointType qp : SurfaceQuadratureType( gridPart(), intersection, 2*basis.order(), SurfaceQuadratureType::INSIDE ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.localPosition() );
            values( basis, qp, phi );
            addLinearizedBoundaryIntegrand( integrands, qp, weight, value( u, qp ), phi, basis.size(), j );
          }
        }

        // addSkeletonIntegral

        template< bool conforming, class Integrands, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Integrands &integrands, const Intersection &intersection, const U &uIn, const U &uOut, W &wIn ) const
        {
          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, 2*wIn.order(), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const RangeFieldType weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];

            addSkeletonIntegrand( integrands, qpIn, qpOut, weight, value( uIn, qpIn ), value( uOut, qpOut ), wIn );
          }
        }

        template< bool conforming, class Integrands, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Integrands &integrands, const Intersection &intersection, const U &uIn, const U &uOut, W &wIn, W &wOut ) const
        {
          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, 2*std::max( wIn.order(), wOut.order() ), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const RangeFieldType weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];

            addSkeletonIntegrand( integrands, qpIn, qpOut, weight, value( uIn, qpIn ), value( uOut, qpOut ), wIn, wOut );
          }
        }

        template< class Integrands, class Intersection, class U, class... W >
        void addSkeletonIntegral ( Integrands &integrands, const Intersection &intersection, const U &uIn, const U &uOut, W &... w ) const
        {
          if( !integrands.init( intersection ) )
            return;

          if( intersection.conforming() )
            addSkeletonIntegral< true >( integrands, intersection, uIn, uOut, w... );
          else
            addSkeletonIntegral< false >( integrands, intersection, uIn, uOut, w... );
        }

        template< bool conforming, class Integrands, class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( const Integrands &integrands, const Intersection &intersection, const U &uIn, const U &uOut, ValueVectorType &phiIn, ValueVectorType &phiOut, J &... j ) const
        {
          const auto &basisIn = discreteFunctionSpace().basisFunctionSet( uIn.entity() );
          const auto &basisOut = discreteFunctionSpace().basisFunctionSet( uOut.entity() );

          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, 2*std::max( basisIn.order(), basisOut.order() ), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const RangeFieldType weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];

            values( basisIn, qpIn, phiIn );
            values( basisOut, qpOut, phiOut );

            addLinearizedSkeletonIntegrand( integrands, qpIn, qpOut, weight, value( uIn, qpIn ), value( uOut, qpOut ), phiIn, basisIn.size(), phiOut, basisOut.size(), j... );
          }
        }

        template< class Integrands, class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( Integrands &integrands, const Intersection &intersection, const U &uIn, const U &uOut, ValueVectorType &phiIn, ValueVectorType &phiOut, J &... j ) const
        {
          if( !integrands.init( intersection ) )
            return;

          if( intersection.conforming() )
            addLinearizedSkeletonIntegral< true >( integrands, intersection, uIn, uOut, phiIn, phiOut, j... );
          else
            addLinearizedSkeletonIntegral< false >( integrands, intersection, uIn, uOut, phiIn, phiOut, j... );
        }

        // constructor

        explicit GalerkinOperator ( const DiscreteFunctionSpaceType &dfSpace ) : dfSpace_( dfSpace ) {}

        // evaluate

      private:
        template< class Integrands, class GridFunction, class DiscreteFunction, std::enable_if_t< !HasSkeletonIntegrand< Integrands >::value, int > = 0 >
        void evaluate ( Integrands &integrands, const GridFunction &u, DiscreteFunction &w, std::false_type ) const
        {
          w.clear();

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wLocal( discreteFunctionSpace() );

          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uLocal = u.localFunction( entity );

            wLocal.init( entity );
            wLocal.clear();

            if( hasInteriorIntegrand( integrands ) )
              addInteriorIntegral( integrands, uLocal, wLocal );

            if( hasBoundaryIntegrand( integrands ) && entity.hasBoundaryIntersections() )
            {
              for( const auto &intersection : intersections( gridPart(), entity ) )
              {
                if( intersection.boundary() )
                  addBoundaryIntegral( integrands, intersection, uLocal, wLocal );
              }
            }

            w.addLocalDofs( entity, wLocal.localDofVector() );
          }

          w.communicate();
        }

        template< class Integrands, class GridFunction, class DiscreteFunction >
        void evaluate ( Integrands &integrands, const GridFunction &u, DiscreteFunction &w, std::true_type ) const
        {
          w.clear();

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wInside( discreteFunctionSpace() ), wOutside( discreteFunctionSpace() );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uInside = u.localFunction( inside );

            wInside.init( inside );
            wInside.clear();

            if( hasInteriorIntegrand( integrands ) )
              addInteriorIntegral( integrands, uInside, wInside );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              if( intersection.boundary() )
              {
                if( hasBoundaryIntegrand( integrands ) )
                  addBoundaryIntegral( integrands, intersection, uInside, wInside );
              }
              else if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                if( outside.partitionType() != InteriorEntity )
                  addSkeletonIntegral( integrands, intersection, uInside, u.localFunction( outside ), wInside );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  wOutside.init( outside );
                  wOutside.clear();
                  addSkeletonIntegral( integrands, intersection, uInside, u.localFunction( outside ), wInside, wOutside );
                  w.addLocalDofs( outside, wOutside.localDofVector() );
                }
              }
            }

            w.addLocalDofs( inside, wInside.localDofVector() );
          }
          w.communicate();
        }

      public:
        template< class Integrands, class GridFunction, class DiscreteFunction >
        void evaluate ( Integrands &integrands, const GridFunction &u, DiscreteFunction &w ) const
        {
          if( hasSkeletonIntegrand( integrands ) )
            evaluate( integrands, u, w, std::true_type() );
          else
            evaluate( integrands, u, w, std::false_type() );
        }

        // assemble

      private:
        template< class Integrands, class GridFunction, class JacobianOperator >
        void assemble ( Integrands &integrands, const GridFunction &u, JacobianOperator &jOp, std::false_type ) const
        {
          DiagonalStencil< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > stencil( discreteFunctionSpace(), discreteFunctionSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = discreteFunctionSpace().blockMapper().maxNumDofs() * discreteFunctionSpace().localBlockSize;
          ValueVectorType phi( maxNumLocalDofs, maxNumLocalDofs );

          TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > jOpLocal( discreteFunctionSpace(), discreteFunctionSpace() );

          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uLocal = u.localFunction( entity );

            jOpLocal.init( entity, entity );
            jOpLocal.clear();

            if( hasInteriorIntegrand( integrands ) )
              addLinearizedInteriorIntegral( integrands, uLocal, phi, jOpLocal );

            if( hasBoundaryIntegrand( integrands ) && entity.hasBoundaryIntersections() )
            {
              for( const auto &intersection : intersections( gridPart(), entity ) )
              {
                if( intersection.boundary() )
                  addLinearizedBoundaryIntegral( integrands, intersection, uLocal, phi, jOpLocal );
              }
            }

            jOp.addLocalMatrix( entity, entity, jOpLocal );
          }

          jOp.communicate();
        }

        template< class Integrands, class GridFunction, class JacobianOperator >
        void assemble ( Integrands &integrands, const GridFunction &u, JacobianOperator &jOp, std::true_type ) const
        {
          typedef TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > TemporaryLocalMatrix;

          DiagonalAndNeighborStencil< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > stencil( discreteFunctionSpace(), discreteFunctionSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = discreteFunctionSpace().blockMapper().maxNumDofs() * discreteFunctionSpace().localBlockSize;
          ValueVectorType phiIn( maxNumLocalDofs, maxNumLocalDofs ), phiOut( maxNumLocalDofs, maxNumLocalDofs );

          TemporaryLocalMatrix jOpInIn( discreteFunctionSpace(), discreteFunctionSpace() ), jOpOutIn( discreteFunctionSpace(), discreteFunctionSpace() );
          TemporaryLocalMatrix jOpInOut( discreteFunctionSpace(), discreteFunctionSpace() ), jOpOutOut( discreteFunctionSpace(), discreteFunctionSpace() );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uIn = u.localFunction( inside );

            jOpInIn.init( inside, inside );
            jOpInIn.clear();

            if( hasInteriorIntegrand( integrands ) )
              addLinearizedInteriorIntegral( integrands, uIn, phiIn, jOpInIn );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              if( intersection.boundary() )
              {
                if( hasBoundaryIntegrand( integrands ) )
                  addLinearizedBoundaryIntegral( integrands, intersection, uIn, phiIn, jOpInIn );
              }
              else if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                jOpOutIn.init( outside, inside );
                jOpOutIn.clear();

                if( outside.partitionType() != InteriorEntity )
                  addLinearizedSkeletonIntegral( integrands, intersection, uIn, u.localFunction( outside ), phiIn, phiOut, jOpInIn, jOpOutIn );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  jOpInOut.init( inside, outside );
                  jOpInOut.clear();
                  jOpOutOut.init( outside, outside );
                  jOpOutOut.clear();

                  addLinearizedSkeletonIntegral( integrands, intersection, uIn, u.localFunction( outside ), phiIn, phiOut, jOpInIn, jOpOutIn, jOpInOut, jOpOutOut );

                  jOp.addLocalMatrix( inside, outside, jOpInOut );
                  jOp.addLocalMatrix( outside, outside, jOpOutOut );
                }

                jOp.addLocalMatrix( outside, inside, jOpOutIn );
              }
            }

            jOp.addLocalMatrix( inside, inside, jOpInIn );
          }

          jOp.communicate();
        }

      public:
        template< class Integrands, class GridFunction, class JacobianOperator >
        void assemble ( Integrands &integrands, const GridFunction &u, JacobianOperator &jOp ) const
        {
          if( hasSkeletonIntegrand( integrands ) )
            assemble( integrands, u, jOp, std::true_type() );
          else
            assemble( integrands, u, jOp, std::false_type() );
        }

        // accessors

        const DiscreteFunctionSpaceType &discreteFunctionSpace () const { return dfSpace_; }
        const GridPartType &gridPart () const { return discreteFunctionSpace().gridPart(); }

      private:
        const DiscreteFunctionSpaceType &dfSpace_;
      };

    } // namespace Impl




    // GalerkinOperator
    // ----------------

    template< class DiscreteFunction, class Integrands >
    struct GalerkinOperator
      : public virtual Operator< DiscreteFunction >
    {
      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      template< class... Args >
      GalerkinOperator ( const DiscreteFunctionSpaceType &dfSpace, Args &&... args )
        : impl_( dfSpace ), integrands_( std::forward< Args >( args )... )
      {}

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
      {
        impl_.evaluate( integrands_, u, w );
      }

    protected:
      Impl::GalerkinOperator< DiscreteFunctionSpaceType > impl_;
      mutable Integrands integrands_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class JacobianOperator, class Integrands >
    class DifferentiableGalerkinOperator
      : public GalerkinOperator< typename JacobianOperator::DomainFunctionType, Integrands >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef GalerkinOperator< typename JacobianOperator::DomainFunctionType, Integrands > BaseType;

      template< class, class, SolverType > friend class GalerkinScheme;

    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      template< class... Args >
      DifferentiableGalerkinOperator ( const DiscreteFunctionSpaceType &dfSpace, Args &&... args )
        : BaseType( dfSpace, std::forward< Args >( args )... )
      {}

      void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
      {
        impl_.assemble( integrands_, u, jOp );
      }

    protected:
      using BaseType::impl_;
      using BaseType::integrands_;
    };



    // GalerkinScheme
    // --------------

    template< class Space, class Integrands, SolverType solver >
    struct GalerkinScheme
    {
      typedef Space DiscreteFunctionSpaceType;
      typedef Integrands ModelType;

      typedef Solvers< DiscreteFunctionSpaceType, solver, false > UsedSolverType;
      static_assert( UsedSolverType::solverConfigured, "chosen solver is not configured" );

      typedef typename UsedSolverType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename UsedSolverType::LinearOperatorType LinearOperatorType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      static const int dimRange = DiscreteFunctionSpaceType::dimRange;

      typedef DifferentiableGalerkinOperator< LinearOperatorType, Integrands > GalerkinOperatorType;

      GalerkinScheme ( const DiscreteFunctionSpaceType &dfSpace, const Integrands &integrands, Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container() )
        : galerkinOperator_( dfSpace, integrands ),
          linearOperator_( "assembled Galerkin operator", dfSpace, dfSpace ),
          rhs_( "rhs", dfSpace ),
          parameter_( std::move( parameter ) )
      {}

      void constraint ( const DiscreteFunctionType &u ) {}

      template< class GridFunction >
      void operator() ( const GridFunction &u, DiscreteFunctionType &w )
      {
        galerkinOperator_.impl_.evaluate( galerkinOperator_.integrands_, u, w );
      }

      void solve ( DiscreteFunctionType &solution )
      {
        typedef typename UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
        NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > invOp( galerkinOperator_, parameter_ );
        rhs_.clear();
        invOp( rhs_, solution );
      }

      template< class GridFunction >
      const LinearOperatorType &assemble ( const GridFunction &u )
      {
        galerkinOperator_.impl_.assemble( galerkinOperator_.integrands_, u, *linearOperator_ );
        return *linearOperator_;
      }

      bool mark ( double tolerance ) { return false; }
      double estimate ( const DiscreteFunctionType &solution ) { return 0.0; }

      const GridPartType &gridPart () const { return space().gridPart(); }
      const DiscreteFunctionSpaceType &space() const { return galerkinOperator_.impl_.discreteFunctionSpace(); }

    protected:
      GalerkinOperatorType galerkinOperator_;
      LinearOperatorType linearOperator_;
      mutable DiscreteFunctionType rhs_;
      Dune::Fem::ParameterReader parameter_;
    };



    // DiffusionModelIntegrands
    // ------------------------

    template< class Model >
    struct DiffusionModelIntegrands
    {
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef std::tuple< RangeType, JacobianRangeType > ValueType;

      explicit DiffusionModelIntegrands ( const Model &model ) : model_( &model ) {}

      bool init ( const EntityType &entity ) { return model().init( entity ); }

      bool init ( const IntersectionType &intersection )
      {
        return (intersection.boundary() && model().hasNeumanBoundary() && model().init( intersection.inside() ));
      }

      template< class Point >
      ValueType interior ( const Point &x, const ValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().diffusiveFlux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return std::make_tuple( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linDiffusiveFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return std::make_tuple( source, dFlux );
          };
      }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return std::make_tuple( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return std::make_tuple( alpha, 0 );
          };
      }

      const Model &model () const { return *model_; }

    private:
      const Model *model_;
    };



    // DGDiffusionModelIntegrands
    // --------------------------

    template< class Model >
    struct DGDiffusionModelIntegrands
    {
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef typename FieldTraits< RangeType >::field_type RangeFieldType;

      typedef std::tuple< RangeType, JacobianRangeType > ValueType;

      DGDiffusionModelIntegrands ( const Model &model, RangeFieldType penalty )
        : model_( &model ), penalty_( penalty )
      {}

      bool init ( const EntityType &entity )
      {
        intersection_ = nullptr;
        return model().init( entity );
      }

      bool init ( const IntersectionType &intersection )
      {
        intersection_ = &intersection;
        if( intersection.boundary() )
        {
          const EntityType inside = intersection.inside();
          beta_ = intersection.geometry().volume() / inside.geometry().volume();
          return (model().hasNeumanBoundary() && model().init( inside ));
        }
        else if( intersection.neighbor() )
        {
          const auto volIn = intersection.inside().geometry().volume();
          const auto volOut = intersection.outside().geometry().volume();
          beta_ = intersection.geometry().volume() / std::min( volIn, volOut );
          return true;
        }
      }

      template< class Point >
      ValueType interior ( const Point &x, const ValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().diffusiveFlux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return std::make_tuple( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linDiffusiveFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return std::make_tuple( source, dFlux );
          };
      }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return std::make_tuple( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return std::make_tuple( alpha, 0 );
          };
      }

      template< class Point >
      std::pair< ValueType, ValueType > skeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const EntityType inside = intersection().inside();
        const EntityType outside = intersection().outside();

        const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );
        const auto normal = intersection.unitOuterNormal( xIn.localPosition() );

        ValueType uJump;
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
        {
          std::get< 1 >( uJump )[ i ] = normal;
          std::get< 1 >( uJump )[ i ] *= std::get< 0 >( uJump )[ i ];
        }

        model().init( outside );
        JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uJump ), dFluxPrimeOut );
        model().diffusiveFlux( xOut, 0, 0, dFluxOut );
        dFluxPrimeOut -= dFluxOut;
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uOut ), dFluxOut );

        model().init( inside );
        JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uJump ), dFluxPrimeIn );
        model().diffusiveFlux( xIn, 0, 0, dFluxIn );
        dFluxPrimeIn -= dFluxIn;
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uIn ), dFluxIn );

        RangeType int0 = std::get< 0 >( uJump );
        int0 *= beta_;
        (dFluxIn + dFluxOut).usmv( -half, normal, int0 );

        dFluxPrimeIn *= -half;
        dFluxPrimeOut *= -half;

        return std::make_pair( ValueType( int0, dFluxPrimeOut ), ValueType( -int0, dFluxPrimeOut ) );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const auto normal = intersection.unitOuterNormal( xIn.localPosition() );

        ValueType uJump;
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
        {
          std::get< 1 >( uJump )[ i ] = normal;
          std::get< 1 >( uJump )[ i ] *= std::get< 0 >( uJump )[ i ];
        }

        auto intIn = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiIn ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          ValueType phiJump;
          std::get< 0 >( phiJump ) = std::get< 0 >( phiIn );
          for( int i = 0; i < RangeType::dimension; ++i )
          {
            std::get< 1 >( phiJump )[ i ] = normal;
            std::get< 1 >( phiJump )[ i ] *= -std::get< 0 >( phiJump )[ i ];
          }

          model().init( outside );
          JacobianRangeType dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, 0, std::get< 1 >( phiJump ), dFluxPrimeOut );

          model().init( inside );
          JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiJump ), dFluxPrimeIn );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uIn ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), dFluxIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxIn.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( int0, dFluxPrimeOut ), ValueType( -int0, dFluxPrimeOut ) );
        };

        auto intOut = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiOut ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );
          const auto normal = intersection.unitOuterNormal( xIn.localPosition() );

          ValueType phiJump;
          std::get< 0 >( phiJump ) = std::get< 0 >( phiOut );
          for( int i = 0; i < RangeType::dimension; ++i )
          {
            std::get< 1 >( phiJump )[ i ] = normal;
            std::get< 1 >( phiJump )[ i ] *= std::get< 0 >( phiJump )[ i ];
          }

          model().init( outside );
          JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiJump ), dFluxPrimeOut );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uOut ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), dFluxOut );

          model().init( inside );
          JacobianRangeType dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, 0, std::get< 1 >( phiJump ), dFluxPrimeIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxOut.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( int0, dFluxPrimeOut ), ValueType( -int0, dFluxPrimeOut ) );
        };

        return std::make_pair( intIn, intOut );
      }

      const Model &model () const { return *model_; }

    private:
      const IntersectionType &intersection () const { assert( intersection_ ); return *intersection_; }

      const Model *model_;
      RangeFieldType penalty_;
      const IntersectionType *intersection_ = nullptr;
      RangeFieldType beta_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
