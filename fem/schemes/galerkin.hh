#ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
#define DUNE_FEM_SCHEMES_GALERKIN_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/io/parameter/reader.hh>
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem/schemes/integrands.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      // GalerkinOperator
      // ----------------

      template< class DiscreteFunctionSpace, class Integrands >
      struct GalerkinOperator
      {
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        typedef std::conditional_t< Fem::IntegrandsTraits< Integrands >::isFull, Integrands, FullIntegrands< Integrands > > IntegrandsType;

        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
        typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

        typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

        static_assert( std::is_same< typename Integrands::GridPartType, GridPartType >::value, "DiscreteFunctionSpace and Integrands must be defined on the same grid part." );
        static_assert( std::is_same< typename Integrands::ValueType, std::tuple< RangeType, JacobianRangeType > >::value, "For now, Integrands::ValueType must be std::tuple< RangeType, JacobianRangeType >." );

      private:
        typedef CachingQuadrature< GridPartType, 0 > InteriorQuadratureType;
        typedef QuadraturePointWrapper< InteriorQuadratureType > InteriorQuadraturePointType;

        typedef CachingQuadrature< GridPartType, 1 > SurfaceQuadratureType;
        typedef QuadraturePointWrapper< SurfaceQuadratureType > SurfaceQuadraturePointType;

        template< class... T >
        static std::tuple< std::vector< T >... > makeVectorValue ( std::tuple< T... > );

        typedef typename Integrands::ValueType ValueType;
        typedef decltype( makeVectorValue( std::declval< ValueType >() ) ) ValueVectorType;


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

        template< class QP, class J >
        void addLinearizedSkeletonIntegrand ( const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn ) const
        {
          auto integrand = integrands_.linearizedSkeleton( qpIn, uIn, qpOut, uOut );
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

        template< class QP, class J >
        void addLinearizedSkeletonIntegrand ( const QP &qpIn, const QP &qpOut, const RangeFieldType &weight, const ValueType &uIn, const ValueType &uOut, const ValueVectorType &phiIn, std::size_t colsIn, const ValueVectorType &phiOut, std::size_t colsOut, J &jInIn, J &jOutIn, J &jInOut, J &jOutOut ) const
        {
          auto integrand = integrands_.linearizedSkeleton( qpIn, uIn, qpOut, uOut );
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

      public:
        // interior integral

        template< class U, class W >
        void addInteriorIntegral ( const U &u, W &w ) const
        {
          if( !integrands_.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();
          for( const InteriorQuadraturePointType qp : InteriorQuadratureType( u.entity(), 2*w.order() ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.position() );

            ValueType integrand = integrands_.interior( qp, value( u, qp ) );

            std::get< 0 >( integrand ) *= weight;
            std::get< 1 >( integrand ) *= weight;
            w.axpy( qp, std::get< 0 >( integrand ), std::get< 1 >( integrand ) );
          }
        }

        template< class U, class J >
        void addLinearizedInteriorIntegral ( const U &u, ValueVectorType &phi, J &j ) const
        {
          if( !integrands_.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();
          const auto &basis = discreteFunctionSpace().basisFunctionSet( u.entity() );
          for( const auto qp : InteriorQuadratureType( u.entity(), 2*basis.order() ) )
          {
            const auto weight = qp.weight() * geometry.integrationElement( qp.position() );

            values( basis, qp, phi );
            auto integrand = integrands_.linearizedInterior( qp, value( u, qp ) );

            for( std::size_t col = 0, cols = basis.size(); col < cols; ++col )
            {
              ValueType intPhi = integrand( std::make_tuple( std::get< 0 >( phi )[ col ], std::get< 1 >( phi )[ col ] ) );
              j.column( col ).axpy( std::get< 0 >( phi ), std::get< 1 >( phi ), std::get< 0 >( intPhi ), std::get< 1 >( intPhi ), weight );
            }
          }
        }

        // boundary integral

        template< class Intersection, class U, class W >
        void addBoundaryIntegral ( const Intersection &intersection, const U &u, W &w ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          for( const SurfaceQuadraturePointType qp : SurfaceQuadratureType( gridPart(), intersection, 2*w.order(), SurfaceQuadratureType::INSIDE ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            ValueType integrand = integrands_.boundary( qp, value( u, qp ) );

            std::get< 0 >( integrand ) *= weight;
            std::get< 1 >( integrand ) *= weight;
            w.axpy( qp, std::get< 0 >( integrand ), std::get< 1 >( integrand ) );
          }
        }

        template< class Intersection, class U, class J >
        void addLinearizedBoundaryIntegral ( const Intersection &intersection, const U &u, ValueVectorType &phi, J &j ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          const auto &basis = discreteFunctionSpace().basisFunctionSet( u.entity() );
          for( const SurfaceQuadraturePointType qp : SurfaceQuadratureType( gridPart(), intersection, 2*basis.order(), SurfaceQuadratureType::INSIDE ) )
          {
            const RangeFieldType weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            values( basis, qp, phi );
            auto integrand = integrands_.linearizedBoundary( qp, value( u, qp ) );

            for( std::size_t col = 0, cols = basis.size(); col < cols; ++col )
            {
              ValueType intPhi = integrand( std::make_tuple( std::get< 0 >( phi )[ col ], std::get< 1 >( phi )[ col ] ) );
              j.column( col ).axpy( std::get< 0 >( phi ), std::get< 1 >( phi ), std::get< 0 >( intPhi ), std::get< 1 >( intPhi ), weight );
            }
          }
        }

        // addSkeletonIntegral

      private:
        template< bool conforming, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &wIn ) const
        {
          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, 2*wIn.order(), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const RangeFieldType weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];
            std::pair< ValueType, ValueType > integrand = integrands_.skeleton( qpIn, value( uIn, qpIn ), qpOut, value( uOut, qpOut ) );

            std::get< 0 >( integrand.first ) *= weight;
            std::get< 1 >( integrand.first ) *= weight;
            wIn.axpy( qpIn, std::get< 0 >( integrand.first ), std::get< 1 >( integrand.first ) );
          }
        }

        template< bool conforming, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &wIn, W &wOut ) const
        {
          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, 2*std::max( wIn.order(), wOut.order() ), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const RangeFieldType weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];
            std::pair< ValueType, ValueType > integrand = integrands_.skeleton( qpIn, value( uIn, qpIn ), qpOut, value( uOut, qpOut ) );

            std::get< 0 >( integrand.first ) *= weight;
            std::get< 1 >( integrand.first ) *= weight;
            wIn.axpy( qpIn, std::get< 0 >( integrand.first ), std::get< 1 >( integrand.first ) );

            std::get< 0 >( integrand.second ) *= weight;
            std::get< 1 >( integrand.second ) *= weight;
            wOut.axpy( qpOut, std::get< 0 >( integrand.second ), std::get< 1 >( integrand.second ) );
          }
        }

        template< bool conforming, class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, ValueVectorType &phiIn, ValueVectorType &phiOut, J &... j ) const
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

            addLinearizedSkeletonIntegrand( qpIn, qpOut, weight, value( uIn, qpIn ), value( uOut, qpOut ), phiIn, basisIn.size(), phiOut, basisOut.size(), j... );
          }
        }

      public:
        template< class Intersection, class U, class... W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &... w ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          if( intersection.conforming() )
            addSkeletonIntegral< true >( intersection, uIn, uOut, w... );
          else
            addSkeletonIntegral< false >( intersection, uIn, uOut, w... );
        }

        template< class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, ValueVectorType &phiIn, ValueVectorType &phiOut, J &... j ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          if( intersection.conforming() )
            addLinearizedSkeletonIntegral< true >( intersection, uIn, uOut, phiIn, phiOut, j... );
          else
            addLinearizedSkeletonIntegral< false >( intersection, uIn, uOut, phiIn, phiOut, j... );
        }

        // constructor

        template< class... Args >
        explicit GalerkinOperator ( const DiscreteFunctionSpaceType &dfSpace, Args &&... args )
          : dfSpace_( dfSpace ), integrands_( std::forward< Args >( args )... )
        {}

        // evaluate

      private:
        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, std::false_type ) const
        {
          w.clear();

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wLocal( discreteFunctionSpace() );

          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uLocal = u.localFunction( entity );

            wLocal.init( entity );
            wLocal.clear();

            if( integrands_.hasInterior() )
              addInteriorIntegral( uLocal, wLocal );

            if( integrands_.hasBoundary() && entity.hasBoundaryIntersections() )
            {
              for( const auto &intersection : intersections( gridPart(), entity ) )
              {
                if( intersection.boundary() )
                  addBoundaryIntegral( intersection, uLocal, wLocal );
              }
            }

            w.addLocalDofs( entity, wLocal.localDofVector() );
          }

          w.communicate();
        }

        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, std::true_type ) const
        {
          w.clear();

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wInside( discreteFunctionSpace() ), wOutside( discreteFunctionSpace() );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uInside = u.localFunction( inside );

            wInside.init( inside );
            wInside.clear();

            if( integrands_.hasInterior() )
              addInteriorIntegral( uInside, wInside );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              if( intersection.boundary() )
              {
                if( integrands_.hasBoundary() )
                  addBoundaryIntegral( intersection, uInside, wInside );
              }
              else if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                if( outside.partitionType() != InteriorEntity )
                  addSkeletonIntegral( intersection, uInside, u.localFunction( outside ), wInside );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  wOutside.init( outside );
                  wOutside.clear();
                  addSkeletonIntegral( intersection, uInside, u.localFunction( outside ), wInside, wOutside );
                  w.addLocalDofs( outside, wOutside.localDofVector() );
                }
              }
            }

            w.addLocalDofs( inside, wInside.localDofVector() );
          }

          w.communicate();
        }

      public:
        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w ) const
        {
          if( integrands_.hasSkeleton() )
            evaluate( u, w, std::true_type() );
          else
            evaluate( u, w, std::false_type() );
        }

        // assemble

      private:
        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp, std::false_type ) const
        {
          DiagonalStencil< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > stencil( discreteFunctionSpace(), discreteFunctionSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = discreteFunctionSpace().blockMapper().maxNumDofs() * discreteFunctionSpace().localBlockSize;
          ValueVectorType phi{ std::vector< RangeType >( maxNumLocalDofs ), std::vector< JacobianRangeType >( maxNumLocalDofs ) };

          TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > jOpLocal( discreteFunctionSpace(), discreteFunctionSpace() );

          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uLocal = u.localFunction( entity );

            jOpLocal.init( entity, entity );
            jOpLocal.clear();

            if( integrands_.hasInterior() )
              addLinearizedInteriorIntegral( uLocal, phi, jOpLocal );

            if( integrands_.hasBoundary() && entity.hasBoundaryIntersections() )
            {
              for( const auto &intersection : intersections( gridPart(), entity ) )
              {
                if( intersection.boundary() )
                  addLinearizedBoundaryIntegral( intersection, uLocal, phi, jOpLocal );
              }
            }

            jOp.addLocalMatrix( entity, entity, jOpLocal );
          }

          jOp.communicate();
        }

        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp, std::true_type ) const
        {
          typedef TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > TemporaryLocalMatrix;

          DiagonalAndNeighborStencil< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > stencil( discreteFunctionSpace(), discreteFunctionSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = discreteFunctionSpace().blockMapper().maxNumDofs() * discreteFunctionSpace().localBlockSize;
          ValueVectorType phiIn{ std::vector< RangeType >( maxNumLocalDofs ), std::vector< JacobianRangeType >( maxNumLocalDofs ) };
          ValueVectorType phiOut{ std::vector< RangeType >( maxNumLocalDofs ), std::vector< JacobianRangeType >( maxNumLocalDofs ) };

          TemporaryLocalMatrix jOpInIn( discreteFunctionSpace(), discreteFunctionSpace() ), jOpOutIn( discreteFunctionSpace(), discreteFunctionSpace() );
          TemporaryLocalMatrix jOpInOut( discreteFunctionSpace(), discreteFunctionSpace() ), jOpOutOut( discreteFunctionSpace(), discreteFunctionSpace() );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            const auto uIn = u.localFunction( inside );

            jOpInIn.init( inside, inside );
            jOpInIn.clear();

            if( integrands_.hasInterior() )
              addLinearizedInteriorIntegral( uIn, phiIn, jOpInIn );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              if( intersection.boundary() )
              {
                if( integrands_.hasBoundary() )
                  addLinearizedBoundaryIntegral( intersection, uIn, phiIn, jOpInIn );
              }
              else if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                jOpOutIn.init( outside, inside );
                jOpOutIn.clear();

                if( outside.partitionType() != InteriorEntity )
                  addLinearizedSkeletonIntegral( intersection, uIn, u.localFunction( outside ), phiIn, phiOut, jOpInIn, jOpOutIn );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  jOpInOut.init( inside, outside );
                  jOpInOut.clear();
                  jOpOutOut.init( outside, outside );
                  jOpOutOut.clear();

                  addLinearizedSkeletonIntegral( intersection, uIn, u.localFunction( outside ), phiIn, phiOut, jOpInIn, jOpOutIn, jOpInOut, jOpOutOut );

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
        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp ) const
        {
          if( integrands_.hasSkeleton() )
            assemble( u, jOp, std::true_type() );
          else
            assemble( u, jOp, std::false_type() );
        }

        // accessors

        const DiscreteFunctionSpaceType &discreteFunctionSpace () const { return dfSpace_; }
        const GridPartType &gridPart () const { return discreteFunctionSpace().gridPart(); }

      private:
        const DiscreteFunctionSpaceType &dfSpace_;
        mutable IntegrandsType integrands_;
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
        : impl_( dfSpace, std::forward< Args >( args )... )
      {}

      void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
      {
        impl_.evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, DiscreteFunctionType &w ) const
      {
        return impl_.evaluate( u, w );
      }

      const DiscreteFunctionSpaceType &discreteFunctionSpace () const { return impl_.discreteFunctionSpace(); }

    protected:
      Impl::GalerkinOperator< DiscreteFunctionSpaceType, Integrands > impl_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class JacobianOperator, class Integrands >
    class DifferentiableGalerkinOperator
      : public GalerkinOperator< typename JacobianOperator::DomainFunctionType, Integrands >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef GalerkinOperator< typename JacobianOperator::DomainFunctionType, Integrands > BaseType;

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
        impl_.assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        impl_.assemble( u, jOp );
      }

    protected:
      using BaseType::impl_;
    };



    // AutomaticDifferenceGalerkinOperator
    // -----------------------------------

    template< class DiscreteFunction, class Integrands >
    class AutomaticDifferenceGalerkinOperator
      : public GalerkinOperator< DiscreteFunction, Integrands >,
        public AutomaticDifferenceOperator< DiscreteFunction >
    {
      typedef GalerkinOperator< DiscreteFunction, Integrands > BaseType;
      typedef AutomaticDifferenceOperator< DiscreteFunction > AutomaticDifferenceOperatorType;

    public:
      typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      template< class... Args >
      AutomaticDifferenceGalerkinOperator ( const DiscreteFunctionSpaceType &dfSpace, Args &&... args )
        : BaseType( dfSpace, std::forward< Args >( args )... ), AutomaticDifferenceOperatorType()
      {}
    };



    // ModelDifferentiableGalerkinOperator
    // -----------------------------------

    template < class LinearOperator, class ModelIntegrands >
    struct ModelDifferentiableGalerkinOperator
      : public DifferentiableGalerkinOperator< LinearOperator, ModelIntegrands >
    {
      typedef DifferentiableGalerkinOperator< LinearOperator, ModelIntegrands > BaseType;

      typedef typename ModelIntegrands::ModelType ModelType;

      typedef typename LinearOperator::DomainFunctionType RangeDiscreteFunctionType;
      typedef typename LinearOperator::RangeSpaceType DiscreteFunctionSpaceType;

      ModelDifferentiableGalerkinOperator ( const ModelType &model, const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( dfSpace, model )
      {}

      template< class GridFunction >
      void apply ( const GridFunction &u, RangeDiscreteFunctionType &w ) const
      {
        (*this)( u, w );
      }

      template< class GridFunction >
      void apply ( const GridFunction &u, LinearOperator &jOp ) const
      {
        (*this).jacobian( u, jOp );
      }

      void prepare( RangeDiscreteFunctionType &u ) const {}

      void prepare( const RangeDiscreteFunctionType &u, RangeDiscreteFunctionType &w ) const {}
    };



    // GalerkinScheme
    // --------------

    template< class Integrands, class LinearOperator, class InverseOperator >
    struct GalerkinScheme
    {
      typedef DifferentiableGalerkinOperator< LinearOperator, Integrands > DifferentiableOperatorType;

      typedef typename DifferentiableOperatorType::RangeDiscreteFunctionType DiscreteFunctionType;
      typedef typename DifferentiableOperatorType::JacobianOperatorType LinearOperatorType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      struct SolverInfo
      {
        SolverInfo ( bool converged, int linearIterations, int nonlinearIterations )
          : converged( converged ), linearIterations( linearIterations ), nonlinearIterations( nonlinearIterations )
        {}

        bool converged;
        int linearIterations, nonlinearIterations;
      };

      GalerkinScheme ( const DiscreteFunctionSpaceType &dfSpace, Integrands integrands, ParameterReader parameter = Parameter::container )
        : fullOperator_( dfSpace, std::move( integrands ) ), parameter_( std::move( parameter ) ), linearOperator_( "assembled elliptic operator", dfSpace, dfSpace )
      {}

      const DifferentiableOperatorType &fullOperator() const { return fullOperator_; }

      void constraint ( DiscreteFunctionType &u ) const {}

      template< class GridFunction >
      void operator() ( const GridFunction &u, DiscreteFunctionType &w ) const
      {
        fullOperator()( u, w );
      }

      SolverInfo solve ( DiscreteFunctionType &solution ) const
      {
        DiscreteFunctionType bnd( solution );
        bnd.clear();

        Dune::Fem::NewtonInverseOperator< LinearOperatorType, InverseOperator > invOp( fullOperator(), parameter_ );
        invOp( bnd, solution );

        return SolverInfo( invOp.converged(), invOp.linearIterations(), invOp.iterations() );
      }

      template< class GridFunction >
      const LinearOperatorType &assemble ( const GridFunction &ubar )
      {
        fullOperator().jacobian( ubar, linearOperator_ );
        return linearOperator_;
      }

      bool mark ( double tolerance ) { return false; }
      double estimate ( const DiscreteFunctionType &solution ) { return 0.0; }

      const DiscreteFunctionSpaceType &space () const { return fullOperator_.discreteFunctionSpace(); }
      const GridPartType &gridPart () const { return space().gridPart(); }

    protected:
      DifferentiableOperatorType fullOperator_;
      ParameterReader parameter_;
      LinearOperatorType linearOperator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
