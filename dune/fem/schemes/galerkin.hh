#ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
#define DUNE_FEM_SCHEMES_GALERKIN_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <shared_mutex>
#include <vector>
#include <memory>

#include <dune/common/hybridutilities.hh>
#include <dune/common/timer.hh>

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
#include <dune/fem/common/bindguard.hh>

#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

#include <dune/fem/operator/common/localmatrixcolumn.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/schemes/integrands.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/femscheme.hh>

#include <dune/fem/space/common/capabilities.hh>

// fempy includes
#include <dune/fempy/quadrature/fempyquadratures.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {
      template <class M>
      class CallOrder
      {

        template <class F>
        static int callOrder(const F& f, char)
        {
#ifndef NDEBUG
          std::cerr << "WARNING: no order method available on " << typeid(F).name() << ", defaulting to 1!" << std::endl;
#endif
          return 1;
        }

        template <class F>
        static auto callOrder(const F& f, int) -> decltype( f.order() )
        {
          return f.order();
        }

      public:
        template <class F>
        static int order (const F& f ) { return callOrder(f, 0); }
      };

      // GalerkinOperator
      // ----------------

      template <class Space>
      struct DefaultGalerkinOperatorQuadratureSelector
      {
        typedef typename Space :: GridPartType GridPartType;
        typedef CachingQuadrature< GridPartType, 0, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > InteriorQuadratureType;
        typedef CachingQuadrature< GridPartType, 1, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > SurfaceQuadratureType;
        // typedef CachingQuadrature< GridPartType, 0, Dune::FemPy::FempyQuadratureTraits > InteriorQuadratureType;
        // typedef CachingQuadrature< GridPartType, 1, Dune::FemPy::FempyQuadratureTraits > SurfaceQuadratureType;
      };

      // LocalGalerkinOperator
      // ---------------------

      template< class Integrands, template <class> class QuadSelector = DefaultGalerkinOperatorQuadratureSelector >
      struct LocalGalerkinOperator
      {
        typedef LocalGalerkinOperator<Integrands> ThisType;
        typedef std::conditional_t< Fem::IntegrandsTraits< Integrands >::isFull, Integrands, FullIntegrands< Integrands > > IntegrandsType;

        typedef typename IntegrandsType::GridPartType GridPartType;

        typedef typename GridPartType::ctype ctype;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        // typedef QuadratureSelector
        template <class Space>
        using QuadratureSelector = QuadSelector< Space >;

        // constructor
        template< class... Args >
        explicit LocalGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
          : gridPart_( gridPart ),
            integrands_( std::forward< Args >( args )... ),
            defaultInteriorOrder_( [] (const int order) { return 2 * order; } ),
            defaultSurfaceOrder_ ( [] (const int order) { return 2 * order + 1; } ),
            interiorQuadOrder_(0), surfaceQuadOrder_(0)
        {
        }

      protected:
        typedef typename IntegrandsType::DomainValueType DomainValueType;
        typedef typename IntegrandsType::RangeValueType RangeValueType;
        typedef std::make_index_sequence< std::tuple_size< DomainValueType >::value > DomainValueIndices;
        typedef std::make_index_sequence< std::tuple_size< RangeValueType >::value > RangeValueIndices;


        template< std::size_t... i >
        static auto makeDomainValueVector ( std::size_t maxNumLocalDofs, std::index_sequence< i... > )
        {
          return std::make_tuple( std::vector< std::tuple_element_t< i, DomainValueType > >( maxNumLocalDofs )... );
        }

        static auto makeDomainValueVector ( std::size_t maxNumLocalDofs )
        {
          return makeDomainValueVector( maxNumLocalDofs, DomainValueIndices() );
        }

        template< std::size_t... i >
        static auto makeRangeValueVector ( std::size_t maxNumLocalDofs, std::index_sequence< i... > )
        {
          return std::make_tuple( std::vector< std::tuple_element_t< i, RangeValueType > >( maxNumLocalDofs )... );
        }

        static auto makeRangeValueVector ( std::size_t maxNumLocalDofs )
        {
          return makeRangeValueVector( maxNumLocalDofs, RangeValueIndices() );
        }

        typedef decltype( makeDomainValueVector( 0u ) ) DomainValueVectorType;
        typedef decltype( makeRangeValueVector( 0u ) )  RangeValueVectorType;

        static void resizeDomainValueVector ( DomainValueVectorType& vec, const std::size_t size )
        {
          Hybrid::forEach( DomainValueIndices(), [ &vec, &size ] ( auto i ) {
              std::get< i >( vec ).resize( size );
            } );
        }

        static void resizeRangeValueVector ( RangeValueVectorType& vec, const std::size_t size )
        {
          Hybrid::forEach( RangeValueIndices(), [ &vec, &size ] ( auto i ) {
              std::get< i >( vec ).resize( size );
            } );
        }

      public:
        void prepare( const std::size_t size ) const
        {
          resizeDomainValueVector( phiIn_, size );
          resizeDomainValueVector( phiOut_, size );
          resizeDomainValueVector( basisValues_, size );
          resizeDomainValueVector( domainValues_, size );
        }

        template< class LocalFunction, class Quadrature >
        static void evaluateQuadrature ( const LocalFunction &u, const Quadrature &quad, std::vector< typename LocalFunction::RangeType > &phi )
        {
          u.evaluateQuadrature( quad, phi );
        }

        template< class LocalFunction, class Quadrature>
        static void evaluateQuadrature ( const LocalFunction &u, const Quadrature &quad, std::vector< typename LocalFunction::JacobianRangeType > &phi )
        {
          u.jacobianQuadrature( quad, phi );
        }

        template< class LocalFunction, class Quadrature >
        static void evaluateQuadrature ( const LocalFunction &u, const Quadrature &quad, std::vector< typename LocalFunction::HessianRangeType > &phi )
        {
          u.hessianQuadrature( quad, phi );
        }

      protected:
        template< class LocalFunction, class Point >
        static void value ( const LocalFunction &u, const Point &x, typename LocalFunction::RangeType &phi )
        {
          u.evaluate( x, phi );
        }

        template< class LocalFunction, class Point >
        static void value ( const LocalFunction &u, const Point &x, typename LocalFunction::JacobianRangeType &phi )
        {
          u.jacobian( x, phi );
        }

        template< class LocalFunction, class Point >
        static void value ( const LocalFunction &u, const Point &x, typename LocalFunction::HessianRangeType &phi )
        {
          u.hessian( x, phi );
        }

        template< class LocalFunction, class Point, class... T >
        static void value ( const LocalFunction &u, const Point &x, std::tuple< T... > &phi )
        {
          Hybrid::forEach( std::index_sequence_for< T... >(), [ &u, &x, &phi ] ( auto i ) { LocalGalerkinOperator::value( u, x, std::get< i >( phi ) ); } );
        }

        template< class Basis, class Point >
        static void values ( const Basis &basis, const Point &x, std::vector< typename Basis::RangeType > &phi )
        {
          basis.evaluateAll( x, phi );
        }

        template< class Basis, class Point >
        static void values ( const Basis &basis, const Point &x, std::vector< typename Basis::JacobianRangeType > &phi )
        {
          basis.jacobianAll( x, phi );
        }

        template< class Basis, class Point >
        static void values ( const Basis &basis, const Point &x, std::vector< typename Basis::HessianRangeType > &phi )
        {
          basis.hessianAll( x, phi );
        }

        template< class Basis, class Point, class... T >
        static void values ( const Basis &basis, const Point &x, std::tuple< std::vector< T >... > &phi )
        {
          Hybrid::forEach( std::index_sequence_for< T... >(), [ &basis, &x, &phi ] ( auto i ) { LocalGalerkinOperator::values( basis, x, std::get< i >( phi ) ); } );
        }

        template< class LocalFunction, class Point >
        static DomainValueType domainValue ( const LocalFunction &u, const Point &x )
        {
          DomainValueType phi;
          value( u, x, phi );
          return phi;
        }

        static DomainValueType domainValue ( const unsigned int qpIdx, DomainValueVectorType& vec)
        {
          DomainValueType phi;
          Hybrid::forEach( DomainValueIndices(), [ &qpIdx, &vec, &phi ] ( auto i ) {
              std::get< i > ( phi )  = std::get< i >( vec )[ qpIdx ];
                } );
          return phi;
        }

        template< class LocalFunction, class Quadrature >
        static void domainValue ( const LocalFunction &u, const Quadrature& quadrature, DomainValueVectorType &result  )
        {
          Hybrid::forEach( DomainValueIndices(), [ &u, &quadrature, &result ] ( auto i ) {
              auto& vec = std::get< i >( result );
              vec.resize( quadrature.nop() );
              ThisType::evaluateQuadrature( u, quadrature, vec );
            } );
        }

        template< class Phi, std::size_t... i >
        static auto value ( const Phi &phi, std::size_t col, std::index_sequence< i... > )
        {
          return std::make_tuple( std::get< i >( phi )[ col ]... );
        }

        template< class... T >
        static auto value ( const std::tuple< std::vector< T >... > &phi, std::size_t col )
        {
          return value( phi, col, std::index_sequence_for< T... >() );
        }

        static void assignRange( RangeValueVectorType& ranges, const std::size_t idx, const RangeValueType& range )
        {
          Hybrid::forEach( RangeValueIndices(), [ &ranges, &idx, &range ] ( auto i ) {
              std::get< i >( ranges )[ idx ] = std::get< i >( range );
            });
        }
        template <class W>
        static void assignRange( RangeValueVectorType& ranges, const std::size_t idx, const RangeValueType& range, const W &weight )
        {
          Hybrid::forEach( RangeValueIndices(), [ &ranges, &idx, &range, &weight ] ( auto i ) {
              std::get< i >( ranges )[ idx ]  = std::get< i >( range );
              std::get< i >( ranges )[ idx ] *= weight;
            });
        }

        static void assignDomain( DomainValueVectorType& domains, const std::size_t idx, const DomainValueType& domain )
        {
          Hybrid::forEach( DomainValueIndices(), [ &domains, &idx, &domain ] ( auto i ) {
              std::get< i >( domains )[ idx ] = std::get< i >( domain );
            });
        }

        template <class W, class Quadrature>
        static void axpyQuadrature( W& w, const Quadrature& quadrature, RangeValueVectorType& ranges )
        {
          Hybrid::forEach( RangeValueIndices(), [ &w, &quadrature, &ranges ] ( auto i ) {
              w.axpyQuadrature( quadrature, std::get< i >( ranges ) );
            } );
        }

      public:
        // interior integral

        template< class U, class W >
        void addInteriorIntegral ( const U &u, W &w ) const
        {
          if( !integrands().init( u.entity() ) )
            return;

          const auto& geometry = u.geometry();

          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: InteriorQuadratureType  InteriorQuadratureType;
          const InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder(maxOrder(u, w)) );

          // evaluate u for all quadrature points
          DomainValueVectorType& domains = domainValues_;
          domainValue( u, quadrature, domains );

          auto& ranges = values_;
          resizeRangeValueVector( ranges, quadrature.nop() );

          // evaluate integrands for all quadrature points
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.position() );
            assignRange( ranges, qp.index(), integrands().interior( qp, domainValue( qp.index(), domains ) ), weight );
          }

          // add to w for all quadrature points
          axpyQuadrature( w, quadrature, ranges );
          integrands().unbind();
        }

        template< class U, class J >
        void addLinearizedInteriorIntegral ( const U &u, J &j ) const
        {
          if( !integrands().init( u.entity() ) )
            return;

          const auto &geometry = u.geometry();
          const auto &domainBasis = j.domainBasisFunctionSet();
          const auto &rangeBasis = j.rangeBasisFunctionSet();

          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: InteriorQuadratureType  InteriorQuadratureType;
          const InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder( maxOrder( u, domainBasis, rangeBasis )) );
          const size_t domainSize = domainBasis.size();
          const size_t quadNop = quadrature.nop();

          auto& basisValues = basisValues_;
          resizeDomainValueVector( basisValues, domainSize );

          // evaluate u for all quadrature points
          auto& rangeValues = rangeValues_;
          DomainValueVectorType& domains = domainValues_;
          domainValue( u, quadrature, domains );

          rangeValues.resize( domainSize );
          for( std::size_t col = 0; col < domainSize; ++col )
          {
            resizeRangeValueVector( rangeValues[ col ], quadNop );
          }

          // evaluate all basis functions and integrands
          for( const auto qp : quadrature )
          {
            values( domainBasis, qp, basisValues );
            const auto weight = qp.weight() * geometry.integrationElement( qp.position() );
            auto integrand = integrands().linearizedInterior( qp, domainValue( qp.index(), domains ) );
            for( std::size_t col = 0; col < domainSize; ++col )
            {
              assignRange( rangeValues[ col ], qp.index(), integrand( value( basisValues, col ) ), weight );
            }
          }

          // add to local matrix for all quadrature points and basis functions
          for( std::size_t col = 0; col < domainSize; ++col )
          {
            LocalMatrixColumn< J > jCol( j, col );
            axpyQuadrature( jCol, quadrature, rangeValues[ col ] );
          }
          integrands().unbind();
        }

        // boundary integral

        template< class Intersection, class U, class W >
        void addBoundaryIntegral ( const Intersection &intersection, const U &u, W &w ) const
        {
          if( !integrands().init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          const SurfaceQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder( u, w )), SurfaceQuadratureType::INSIDE );
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            RangeValueType integrand = integrands().boundary( qp, domainValue( u, qp ) );

            Hybrid::forEach( RangeValueIndices(), [ &qp, &w, &integrand, weight ] ( auto i ) {
                std::get< i >( integrand ) *= weight;
                w.axpy( qp, std::get< i >( integrand ) );
              } );
          }
          integrands().unbind();
        }

        template< class Intersection, class U, class J >
        void addLinearizedBoundaryIntegral ( const Intersection &intersection, const U &u, J &j ) const
        {
          if( !integrands().init( intersection ) )
            return;

          DomainValueVectorType &phi = phiIn_;

          const auto geometry = intersection.geometry();
          const auto &domainBasis = j.domainBasisFunctionSet();
          const auto &rangeBasis = j.rangeBasisFunctionSet();

          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          const SurfaceQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder(u, domainBasis, rangeBasis )), SurfaceQuadratureType::INSIDE );
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            values( domainBasis, qp, phi );
            auto integrand = integrands().linearizedBoundary( qp, domainValue( u, qp ) );

            for( std::size_t col = 0, cols = domainBasis.size(); col < cols; ++col )
            {
              LocalMatrixColumn< J > jCol( j, col );
              RangeValueType intPhi = integrand( value( phi, col ) );

              Hybrid::forEach( RangeValueIndices(), [ &qp, &jCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi ) *= weight;
                  jCol.axpy( qp, std::get< i >( intPhi ) );
                } );
            }
          }
          integrands().unbind();
        }

        // addSkeletonIntegral

      protected:
        template< bool conforming, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &wIn ) const
        {
          const auto geometry = intersection.geometry();

          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          typedef IntersectionQuadrature< SurfaceQuadratureType, conforming > IntersectionQuadratureType;
          const IntersectionQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder( uIn, uOut, wIn)), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const ctype weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];
            std::pair< RangeValueType, RangeValueType > integrand = integrands().skeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );

            Hybrid::forEach( RangeValueIndices(), [ &qpIn, &wIn, &integrand, weight ] ( auto i ) {
                std::get< i >( integrand.first ) *= weight;
                wIn.axpy( qpIn, std::get< i >( integrand.first ) );
              } );
          }
        }

        template< bool conforming, class Intersection, class U, class W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &wIn, W &wOut ) const
        {
          const auto geometry = intersection.geometry();
          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          typedef IntersectionQuadrature< SurfaceQuadratureType, conforming > IntersectionQuadratureType;
          const IntersectionQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder( uIn, uOut, wIn, wOut)), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const ctype weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];
            std::pair< RangeValueType, RangeValueType > integrand = integrands().skeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );

            Hybrid::forEach( RangeValueIndices(), [ &qpIn, &wIn, &qpOut, &wOut, &integrand, weight ] ( auto i ) {
                std::get< i >( integrand.first ) *= weight;
                wIn.axpy( qpIn, std::get< i >( integrand.first ) );

                std::get< i >( integrand.second ) *= weight;
                wOut.axpy( qpOut, std::get< i >( integrand.second ) );
              } );
          }
        }

        template< bool conforming, class Intersection, class U, class J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection,
                                             const U &uIn, const U &uOut, J &jInIn, J &jOutIn ) const
        {
          DomainValueVectorType &phiIn  = phiIn_;
          DomainValueVectorType &phiOut = phiOut_;

          const auto &domainBasisIn = jInIn.domainBasisFunctionSet();
          const auto &domainBasisOut = jOutIn.domainBasisFunctionSet();

          const auto &rangeBasisIn = jInIn.rangeBasisFunctionSet();

          const int order = std::max( maxOrder(uIn, uOut), maxOrder( domainBasisIn, domainBasisOut, rangeBasisIn ));

          const auto geometry = intersection.geometry();
          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          typedef IntersectionQuadrature< SurfaceQuadratureType, conforming > IntersectionQuadratureType;
          const IntersectionQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(order), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const ctype weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];

            values( domainBasisIn, qpIn, phiIn );
            values( domainBasisOut, qpOut, phiOut );

            auto integrand = integrands().linearizedSkeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );
            for( std::size_t col = 0, cols = domainBasisIn.size(); col < cols; ++col )
            {
              LocalMatrixColumn< J > jInInCol( jInIn, col );
              std::pair< RangeValueType, RangeValueType > intPhi = integrand.first( value( phiIn, col ) );

              Hybrid::forEach( RangeValueIndices(), [ &qpIn, &jInInCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jInInCol.axpy( qpIn, std::get< i >( intPhi.first ) );
                } );
            }
            for( std::size_t col = 0, cols = domainBasisOut.size(); col < cols; ++col )
            {
              LocalMatrixColumn< J > jOutInCol( jOutIn, col );
              std::pair< RangeValueType, RangeValueType > intPhi = integrand.second( value( phiOut, col ) );

              Hybrid::forEach( RangeValueIndices(), [ &qpIn, &jOutInCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jOutInCol.axpy( qpIn, std::get< i >( intPhi.first ) );
                } );
            }
          }
        }

        template< bool conforming, class Intersection, class U, class J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut,
                                             J &jInIn, J &jOutIn, J &jInOut, J &jOutOut ) const
        {
          DomainValueVectorType &phiIn  = phiIn_;
          DomainValueVectorType &phiOut = phiOut_;

          const auto &domainBasisIn = jInIn.domainBasisFunctionSet();
          const auto &domainBasisOut = jOutIn.domainBasisFunctionSet();

          const auto &rangeBasisIn = jInIn.rangeBasisFunctionSet();
          const auto &rangeBasisOut = jInOut.rangeBasisFunctionSet();

          const int order = std::max( maxOrder(uIn, uOut), maxOrder( domainBasisIn, domainBasisOut, rangeBasisIn, rangeBasisOut ));

          const auto geometry = intersection.geometry();
          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          typedef IntersectionQuadrature< SurfaceQuadratureType, conforming > IntersectionQuadratureType;
          const IntersectionQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(order), false );
          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const ctype weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature.inside()[ qp ];
            const auto qpOut = quadrature.outside()[ qp ];

            values( domainBasisIn, qpIn, phiIn );
            values( domainBasisOut, qpOut, phiOut );

            auto integrand = integrands().linearizedSkeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );
            for( std::size_t col = 0, cols = domainBasisIn.size(); col < cols; ++col )
            {
              LocalMatrixColumn< J > jInInCol( jInIn, col );
              LocalMatrixColumn< J > jInOutCol( jInOut, col );
              std::pair< RangeValueType, RangeValueType > intPhi = integrand.first( value( phiIn, col ) );

              Hybrid::forEach( RangeValueIndices(), [ &qpIn, &jInInCol, &qpOut, &jInOutCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jInInCol.axpy( qpIn, std::get< i >( intPhi.first ) );

                  std::get< i >( intPhi.second ) *= weight;
                  jInOutCol.axpy( qpOut, std::get< i >( intPhi.second ) );
                } );
            }
            for( std::size_t col = 0, cols = domainBasisOut.size(); col < cols; ++col )
            {
              LocalMatrixColumn< J > jOutInCol( jOutIn, col );
              LocalMatrixColumn< J > jOutOutCol( jOutOut, col );
              std::pair< RangeValueType, RangeValueType > intPhi = integrand.second( value( phiOut, col ) );

              Hybrid::forEach( RangeValueIndices(), [ &qpIn, &jOutInCol, &qpOut, &jOutOutCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jOutInCol.axpy( qpIn, std::get< i >( intPhi.first ) );

                  std::get< i >( intPhi.second ) *= weight;
                  jOutOutCol.axpy( qpOut, std::get< i >( intPhi.second ) );
                } );
            }
          }
        }

      public:
        template< class Intersection, class U, class... W >
        void addSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, W &... w ) const
        {
          if( !integrands().init( intersection ) )
            return;

          if( intersection.conforming() )
            addSkeletonIntegral< true >( intersection, uIn, uOut, w... );
          else
            addSkeletonIntegral< false >( intersection, uIn, uOut, w... );
          integrands().unbind();
        }

        template< class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, J &... j ) const
        {
          if( !integrands().init( intersection ) )
            return;

          if( intersection.conforming() )
            addLinearizedSkeletonIntegral< true >( intersection, uIn, uOut, j... );
          else
            addLinearizedSkeletonIntegral< false >( intersection, uIn, uOut, j... );
          integrands().unbind();
        }

        void setQuadratureOrders(unsigned int interior, unsigned int surface)
        {
          interiorQuadOrder_ = interior;
          surfaceQuadOrder_  = surface;
        }

        void setQuadratureOrderFunctions( std::function<int(const int)> interiorOrder,
                                          std::function<int(const int)> surfaceOrder )
        {
          defaultInteriorOrder_ = interiorOrder;
          defaultSurfaceOrder_  = surfaceOrder;
        }

        IntegrandsType& model() const
        {
          return integrands();
        }
        bool nonlinear()   const { return model().nonlinear(); }
        bool hasInterior() const { return model().hasInterior(); }
        bool hasSkeleton() const { return model().hasSkeleton(); }
        bool hasBoundary() const { return model().hasBoundary(); }

      private:
        IntegrandsType& integrands() const
        {
          return integrands_;
        }

      public:
        // accessors
        const GridPartType &gridPart () const { return gridPart_; }

        unsigned int interiorQuadratureOrder(unsigned int order) const { return interiorQuadOrder_ == 0 ? defaultInteriorOrder_(order) : interiorQuadOrder_; }
        unsigned int surfaceQuadratureOrder(unsigned int order)  const { return surfaceQuadOrder_  == 0 ? defaultSurfaceOrder_ (order) : surfaceQuadOrder_;  }

      protected:
        template <class U>
        int maxOrder( const U& u ) const
        {
          return CallOrder< U > :: order( u );
        }

        template< class U, class W >
        int maxOrder( const U& u, const W& w ) const
        {
          return std::max( maxOrder( u ), maxOrder( w ) );
        }

        template< class U, class V, class W >
        int maxOrder( const U& u, const V& v, const W& w ) const
        {
          return std::max( maxOrder( u, v ), maxOrder( w ) );
        }

        template< class U, class V, class W, class X >
        int maxOrder( const U& u, const V& v, const W& w, const X& x ) const
        {
          return std::max( maxOrder( u, v ), maxOrder( w, x) );
        }

      protected:
        const GridPartType &gridPart_;

        mutable IntegrandsType integrands_;

        mutable std::function<int(const int)> defaultInteriorOrder_;
        mutable std::function<int(const int)> defaultSurfaceOrder_;

        unsigned int interiorQuadOrder_;
        unsigned int surfaceQuadOrder_;

        mutable std::vector< RangeValueVectorType > rangeValues_;
        mutable RangeValueVectorType   values_;
        mutable DomainValueVectorType  phiIn_;
        mutable DomainValueVectorType  phiOut_;
        mutable DomainValueVectorType  basisValues_;
        mutable DomainValueVectorType  domainValues_;
      };



      // GalerkinOperator
      // ----------------

      template< class GridPart >
        // Integrands, template <class> class QuadSelector = DefaultGalerkinOperatorQuadratureSelector >
      struct GalerkinOperator
      {
        typedef GridPart GridPartType;
        typedef GalerkinOperator< GridPartType > ThisType;

        typedef typename GridPartType::ctype ctype;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        // constructor
        explicit GalerkinOperator ( const GridPartType &gridPart )
          : gridPart_( gridPart ),
            gridSizeInterior_( 0 )
        {
        }

      protected:
        template <class IntegrandsTuple>
        bool hasBoundary( const IntegrandsTuple& integrandsTuple ) const
        {
          typedef std::make_index_sequence< std::tuple_size< IntegrandsTuple >::value > Indices;
          bool hasBoundary = false ;
          Hybrid::forEach( Indices(), [&integrandsTuple, &hasBoundary]( auto i ) {
                if( std::get< i > (integrandsTuple).hasBoundary() )
                {
                  hasBoundary = true ;
                  return ;
                }
              });
          return hasBoundary;
        }

        template< class GridFunction, class DiscreteFunction, class Iterators, class IntegrandsTuple, class Functor, bool hasSkeleton >
        void evaluateImpl ( const GridFunction &u, DiscreteFunction &w, const Iterators& iterators,
                                   const IntegrandsTuple& integrandsTuple, Functor& addLocalDofs, std::integral_constant<bool, hasSkeleton> ) const
        {
          Dune::Fem::ConstLocalFunction< GridFunction > uInside( u );
          Dune::Fem::ConstLocalFunction< GridFunction > uOutside( u );

          typedef typename DiscreteFunction::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
          TemporaryLocalFunction< DiscreteFunctionSpaceType > wInside( w.space() ), wOutside( w.space() );

          // element counter
          gridSizeInterior_ = 0;

          typedef std::make_index_sequence< std::tuple_size< IntegrandsTuple >::value > Indices;

          // true if one of the integrands has a boundary term
          const bool hasBnd = hasBoundary( integrandsTuple );

          const auto &indexSet = gridPart().indexSet();
          const auto end = iterators.end();
          for( auto it = iterators.begin(); it != end; ++it )
          {
            // assert( iterators.thread( *it ) == MPIManager::thread() );
            const EntityType inside = *it ;

            // increase counter for interior elements
            ++gridSizeInterior_;

            auto uGuard = bindGuard( uInside, inside );
            auto wGuard = bindGuard( wInside, inside );
            wInside.clear();

            auto addInteriorIntegral = [&integrandsTuple, &uInside, &wInside]( auto i )
            {
              const auto& integrands = std::get< i >( integrandsTuple );
              if( integrands.hasInterior() )
                integrands.addInteriorIntegral( uInside, wInside );
            };
            // add interior integral of any integrands
            Hybrid::forEach( Indices(), addInteriorIntegral );

            if( hasSkeleton || (hasBnd && inside.hasBoundaryIntersections() ) )
            {
              for( const auto &intersection : intersections( gridPart(), inside ) )
              {
                bool neighbor = false;
                if constexpr ( hasSkeleton )
                {
                  // check neighbor first since on periodic boundaries both,
                  // neighbor and boundary are true, so we treat neighbor first
                  if( intersection.neighbor() )
                  {
                    neighbor = true;
                    const EntityType outside = intersection.outside();

                    if( outside.partitionType() != InteriorEntity )
                    {
                      auto uOutGuard = bindGuard( uOutside, outside );

                      auto addSkeletonIntegral = [&integrandsTuple, &intersection, &uInside, &uOutside, &wInside] ( auto i )
                      {
                        const auto& integrands = std::get< i >( integrandsTuple );
                        if( integrands.hasSkeleton() )
                          integrands.addSkeletonIntegral( intersection, uInside, uOutside, wInside );
                      };
                      // add skeleton integral of any integrands
                      Hybrid::forEach( Indices(), addSkeletonIntegral );
                    }
                    else if( indexSet.index( inside ) < indexSet.index( outside ) )
                    {
                      auto uOutGuard = bindGuard( uOutside, outside );
                      auto wOutGuard = bindGuard( wOutside, outside );
                      wOutside.clear();

                      auto addSkeletonIntegral = [&integrandsTuple, &intersection, &uInside, &uOutside, &wInside, &wOutside] ( auto i )
                      {
                        const auto& integrands = std::get< i >( integrandsTuple );
                        if( integrands.hasSkeleton() )
                          integrands.addSkeletonIntegral( intersection, uInside, uOutside, wInside, wOutside );
                      };
                      // add skeleton integral of any integrands
                      Hybrid::forEach( Indices(), addSkeletonIntegral );

                      // addLocalDofs calls w.addLocalDofs but also
                      // prevents race condition for thread parallel runs
                      addLocalDofs( outside, wOutside );
                    }
                  }
                } // end skeleton

                if( ! neighbor && intersection.boundary() )
                {
                  auto addBoundaryIntegral = [&integrandsTuple, &intersection, &uInside, &wInside]( auto i )
                  {
                    const auto& integrands = std::get< i >( integrandsTuple );
                    if( integrands.hasBoundary() )
                      integrands.addBoundaryIntegral( intersection, uInside, wInside );
                  };
                  // add boundary integral of any integrands
                  Hybrid::forEach( Indices(), addBoundaryIntegral );
                } // end boundary
              }
            } // end intersections

            addLocalDofs( inside, wInside );
          }
        }

        template <class Space>
        struct InsideEntity
        {
          typedef typename Space::EntityType EntityType;
          template <class Iterators>
          InsideEntity(const Space &space, const Iterators& iterators)
          : space_(space), dofThread_(space.size(),-1)
          , thread_( iterators.thread() )
          {
            const auto& mapper = space_.blockMapper();
            for (const auto &entity : space_)
            {
              int t=iterators.threadParallel(entity);
              mapper.mapEach(entity, [ this, t ] ( int local, auto global )
                { dofThread_[global] = (dofThread_[global]==t || dofThread_[global]==-1)?
                                       t : -2 ; } ); // -2: shared dof
            }
          }
          bool operator()(const EntityType &entity) const
          {
            bool needsLocking = false;
            space_.blockMapper().mapEach(entity,
                [ this, &needsLocking ] ( int local, auto global )
                { needsLocking = (needsLocking || dofThread_[global]!=thread_); });
            return !needsLocking;
          }
          const Space &space_;
          std::vector<int> dofThread_;
          int thread_;
        };

        template <class DiscreteFunction>
        struct AddLocalEvaluate
        {
          AddLocalEvaluate(DiscreteFunction &w)
          : w_(w) {}
          template <class LocalDofs>
          void operator () (const EntityType& entity, const LocalDofs& wLocal ) const
          {
            w_.addLocalDofs( entity, wLocal.localDofVector() );
          }
          DiscreteFunction &w_;
        };

        template <class DiscreteFunction>
        struct AddLocalEvaluateLocked : public AddLocalEvaluate<DiscreteFunction>
        {
          typedef AddLocalEvaluate<DiscreteFunction> BaseType;

          std::shared_mutex& mutex_;
          InsideEntity<typename DiscreteFunction::DiscreteFunctionSpaceType> inside_;

          template <class Iterators>
          AddLocalEvaluateLocked(DiscreteFunction &w, std::shared_mutex& mtx, const Iterators &iterators)
          : BaseType(w), mutex_(mtx), inside_(w.space(),iterators) {}

          template <class LocalDofs>
          void operator () (const EntityType& entity, const LocalDofs& wLocal ) const
          {
            // call addLocalDofs on w
            if (inside_(entity))
            {
              std::shared_lock<std::shared_mutex> guard ( mutex_ );
              BaseType::operator()( entity, wLocal );
            }
            else
            {
              // lock mutex (unlock on destruction)
              std::lock_guard<std::shared_mutex> guard ( mutex_ );
              BaseType::operator()( entity, wLocal );
            }
          }
        };

        template< class GridFunction, class DiscreteFunction, class Iterators, class IntegrandsTuple, class Functor >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, const Iterators& iterators,
                        const IntegrandsTuple& integrandsTuple, Functor& addLocalDofs ) const
        {
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename DiscreteFunction::GridPartType, GridPartType >::value, "Argument 'w' and Integrands must be defined on the same grid part." );

          if( hasSkeleton( integrandsTuple ) )
            evaluateImpl( u, w, iterators, integrandsTuple, addLocalDofs, std::true_type() );
          else
            evaluateImpl( u, w, iterators, integrandsTuple, addLocalDofs, std::false_type() );
        }

      public:
        template <class IntegrandsTuple>
        bool hasSkeleton( const IntegrandsTuple& integrandsTuple ) const
        {
          typedef std::make_index_sequence< std::tuple_size< IntegrandsTuple >::value > Indices;
          bool hasSkeleton = false ;
          Hybrid::forEach( Indices(), [&integrandsTuple, &hasSkeleton] ( auto i ) {
                            if( std::get< i >( integrandsTuple ).hasSkeleton() )
                            {
                              hasSkeleton = true;
                              return ;
                            }
                         });
          return hasSkeleton ;
        }

        template< class GridFunction, class DiscreteFunction, class Iterators, class IntegrandsTuple >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, const Iterators& iterators,
                        const IntegrandsTuple& integrandsTuple, std::shared_mutex& mtx ) const
        {
          AddLocalEvaluateLocked<DiscreteFunction> addLocalEvaluate(w,mtx,iterators);
          evaluate( u, w, iterators, integrandsTuple, addLocalEvaluate );
        }

        template< class GridFunction, class DiscreteFunction, class Iterators, class IntegrandsTuple >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, const Iterators& iterators, const IntegrandsTuple& integrandsTuple ) const
        {
          AddLocalEvaluate<DiscreteFunction> addLocalEvaluate(w);
          evaluate( u, w, iterators, integrandsTuple, addLocalEvaluate );
        }

      protected:
        template<class T, int length>
        class FiniteStack
        {
        public :
          // Makes empty stack
          FiniteStack () : _f(0) {}

          // Returns true if the stack is empty
          bool empty () const { return _f <= 0; }

          // Returns true if the stack is full
          bool full () const { return (_f >= length); }

          // clear stack
          void clear() { _f = 0; }

          // Puts a new object onto the stack
          void push (const T& t)
          {
            assert ( _f < length );
            _s[_f++] = t;
          }

          // Removes and returns the uppermost object from the stack
          T pop () {
            assert ( _f > 0 );
            return _s[--_f];
          }

          // Returns the uppermost object on the stack
          T top () const {
            assert ( _f > 0 );
            return _s[_f-1];
          }

          // stacksize
          int size () const { return _f; }

        private:
          T   _s[length]; // the stack
          int _f;         // actual position in stack
        };


        template <class JacobianOperator>
        struct AddLocalAssemble
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;
          typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
          JacobianOperator &jOp_;
          std::vector< TemporaryLocalMatrixType > jOpLocal_;

          FiniteStack< TemporaryLocalMatrixType*, 12 > jOpLocalFinalized_;
          FiniteStack< TemporaryLocalMatrixType*, 12 > jOpLocalFree_;

          std::size_t locked, notLocked, timesLocked;
          AddLocalAssemble(JacobianOperator& jOp)
          : jOp_(jOp)
          , jOpLocal_(12, TemporaryLocalMatrixType(jOp_.domainSpace(), jOp_.rangeSpace()))
          , jOpLocalFinalized_()
          , jOpLocalFree_()
          , locked(0), notLocked(0), timesLocked(0)
          {
            for( auto& jOpLocal : jOpLocal_ )
              jOpLocalFree_.push( &jOpLocal );
          }

          TemporaryLocalMatrixType& bind(const EntityType& dE, const EntityType& rE)
          {
            assert( ! jOpLocalFree_.empty() );
            TemporaryLocalMatrixType& lop = *(jOpLocalFree_.pop());
            lop.bind(dE,rE);
            lop.clear();
            return lop;
          }

          void unbind(TemporaryLocalMatrixType &lop)
          {
            notLocked += 1;
            jOp_.addLocalMatrix( lop.domainEntity(), lop.rangeEntity(), lop );
            lop.unbind();
            jOpLocalFree_.push( &lop );
          }

          void finalize()
          {
            locked += jOpLocalFinalized_.size();
            while ( ! jOpLocalFinalized_.empty() )
            {
              TemporaryLocalMatrixType &lop = *(jOpLocalFinalized_.pop());
              jOp_.addLocalMatrix( lop.domainEntity(), lop.rangeEntity(), lop );
              lop.unbind();
              jOpLocalFree_.push( &lop );
            }
          }
        };

        template <class JacobianOperator>
        struct AddLocalAssembleLocked : public AddLocalAssemble<JacobianOperator>
        {
          typedef AddLocalAssemble<JacobianOperator> BaseType;
          typedef typename BaseType::TemporaryLocalMatrixType TemporaryLocalMatrixType;
          using BaseType::jOpLocalFinalized_;
          using BaseType::jOpLocalFree_;

          std::shared_mutex& mutex_;
          InsideEntity<typename JacobianOperator::DomainSpaceType> insideDomain_;
          InsideEntity<typename JacobianOperator::RangeSpaceType>  insideRange_;

          template <class Iterators>
          AddLocalAssembleLocked(JacobianOperator &jOp, std::shared_mutex &mtx, const Iterators &iterators)
          : BaseType(jOp)
          , mutex_(mtx)
          , insideDomain_(jOp.domainSpace(),iterators)
          , insideRange_(jOp.rangeSpace(),iterators)
          {}

          void finalize()
          {
            // lock mutex (unlock on destruction)
            ++BaseType::timesLocked;
            std::lock_guard<std::shared_mutex> guard ( mutex_ );
            BaseType::finalize();
          }

          TemporaryLocalMatrixType& bind(const EntityType& dE, const EntityType& rE)
          {
            if ( jOpLocalFree_.empty() )
            {
              finalize();
            }
            return BaseType::bind(dE,rE);
          }

          void unbind(TemporaryLocalMatrixType &lop)
          {
            /* // always lock
              ++BaseType::timesLocked;
              ++BaseType::locked;
              std::lock_guard guard ( mutex_ );
              BaseType::unbind(lop);
              return;
            */
            if ( insideDomain_(lop.domainEntity()) &&
                 insideRange_(lop.rangeEntity()) )
            {
              std::shared_lock<std::shared_mutex> guard ( mutex_ );
              BaseType::unbind(lop);
            }
            else
            {
              jOpLocalFinalized_.push( &lop );
            }
          }
        };

        template< class GridFunction, class JacobianOperator, class Iterators, class IntegrandsTuple, class Functor, bool hasSkeleton >
        void assembleImpl ( const GridFunction &u, JacobianOperator &jOp, const Iterators& iterators, const IntegrandsTuple& integrandsTuple,
                            Functor& addLocalMatrix, std::integral_constant<bool, hasSkeleton> ) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;

          Dune::Fem::ConstLocalFunction< GridFunction > uIn( u );
          Dune::Fem::ConstLocalFunction< GridFunction > uOut( u );

          typedef std::make_index_sequence< std::tuple_size< IntegrandsTuple >::value > Indices;
          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;

          // initialize local temporary data
          Hybrid::forEach( Indices(), [&integrandsTuple, &maxNumLocalDofs] ( auto i ) {
                           const auto& integrands = std::get< i >( integrandsTuple );
                           integrands.prepare( maxNumLocalDofs );
                         });

          // element counter
          gridSizeInterior_ = 0;

          // true if one of the integrands has a boundary term
          const bool hasBnd = hasBoundary( integrandsTuple );

          const auto &indexSet = gridPart().indexSet();
          // threaded iterators provide from outside
          const auto end = iterators.end();
          for( auto it = iterators.begin(); it != end; ++it )
          {
            // increase counter for interior elements
            ++gridSizeInterior_;

            const EntityType inside = *it;

            auto uiGuard = bindGuard( uIn, inside );

            TemporaryLocalMatrixType& jOpInIn = addLocalMatrix.bind( inside, inside );
            auto addLinearizedInteriorIntegral = [&integrandsTuple, &uIn,  &jOpInIn]( auto i )
            {
              const auto& integrands = std::get< i >( integrandsTuple );
              if( integrands.hasInterior() )
                integrands.addLinearizedInteriorIntegral( uIn, jOpInIn );
            };
            // add interior integral of any integrands
            Hybrid::forEach( Indices(), addLinearizedInteriorIntegral );

            if( hasSkeleton || (hasBnd && inside.hasBoundaryIntersections() ) )
            {
              for( const auto &intersection : intersections( gridPart(), inside ) )
              {
                bool neighbor = false ;
                // check neighbor first since on periodic boundaries both,
                // neighbor and boundary are true, so we treat neighbor first
                if constexpr ( hasSkeleton )
                {
                  if( intersection.neighbor() )
                  {
                    neighbor = true ;
                    const EntityType &outside = intersection.outside();

                    TemporaryLocalMatrixType &jOpOutIn = addLocalMatrix.bind( outside, inside );

                    auto uoGuard = bindGuard( uOut, outside );

                    if( outside.partitionType() != InteriorEntity )
                    {
                      auto addLinearizedSkeletonIntegral = [&integrandsTuple, &intersection, &uIn, &uOut, &jOpInIn, &jOpOutIn]( auto i )
                      {
                        const auto& integrands = std::get< i >( integrandsTuple );
                        if( integrands.hasSkeleton() )
                          integrands.addLinearizedSkeletonIntegral( intersection, uIn, uOut, jOpInIn, jOpOutIn );
                      };
                      // add skeleton integral of any integrands
                      Hybrid::forEach( Indices(), addLinearizedSkeletonIntegral );
                    }
                    else if( indexSet.index( inside ) < indexSet.index( outside ) )
                    {
                      TemporaryLocalMatrixType &jOpInOut = addLocalMatrix.bind( inside, outside );
                      TemporaryLocalMatrixType &jOpOutOut = addLocalMatrix.bind( outside, outside );

                      auto addLinearizedSkeletonIntegral = [&integrandsTuple, &intersection, &uIn, &uOut, &jOpInIn, &jOpOutIn, &jOpInOut, &jOpOutOut]( auto i )
                      {
                        const auto& integrands = std::get< i >( integrandsTuple );
                        if( integrands.hasSkeleton() )
                          integrands.addLinearizedSkeletonIntegral( intersection, uIn, uOut, jOpInIn, jOpOutIn, jOpInOut, jOpOutOut );
                      };
                      // add skeleton integral of any integrands
                      Hybrid::forEach( Indices(), addLinearizedSkeletonIntegral );

                      addLocalMatrix.unbind(jOpInOut);
                      addLocalMatrix.unbind(jOpOutOut);
                    }

                    addLocalMatrix.unbind(jOpOutIn);
                  }
                } // end skeleton

                if( !neighbor && intersection.boundary() )
                {
                  auto addLinearizedBoundaryIntegral = [&integrandsTuple, &intersection, &uIn, &jOpInIn]( auto i )
                  {
                    const auto& integrands = std::get< i >( integrandsTuple );
                    if( integrands.hasBoundary() )
                      integrands.addLinearizedBoundaryIntegral( intersection, uIn, jOpInIn );
                  };
                  // add skeleton integral of any integrands
                  Hybrid::forEach( Indices(), addLinearizedBoundaryIntegral );

                } // end boundary
              }
            } // end intersection
            addLocalMatrix.unbind(jOpInIn);
          }

          // complete the matrix build
          addLocalMatrix.finalize();
        }


        template< class GridFunction, class JacobianOperator, class Iterators, class IntegrandsTuple, class Functor >
        void assemble ( const GridFunction &u, JacobianOperator &jOp, const Iterators& iterators,
                               const IntegrandsTuple& integrandsTuple, Functor& addLocalMatrix, int ) const
        {
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::DomainSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::RangeSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );

          if( hasSkeleton( integrandsTuple ) )
            assembleImpl( u, jOp, iterators, integrandsTuple ,addLocalMatrix, std::true_type() );
          else
            assembleImpl( u, jOp, iterators, integrandsTuple, addLocalMatrix, std::false_type() );
        }

      public:
        template< class GridFunction, class JacobianOperator, class Iterators, class IntegrandsTuple>
        void assemble ( const GridFunction &u, JacobianOperator &jOp, const Iterators& iterators,
                               const IntegrandsTuple& integrandsTuple, std::shared_mutex& mtx) const
        {
          AddLocalAssembleLocked<JacobianOperator> addLocalAssemble( jOp, mtx, iterators);
          assemble( u, jOp, iterators, integrandsTuple, addLocalAssemble, 10 );
          #if 0 // print information about how many times a lock was used during assemble
          std::lock_guard guard ( mtx );
          std::cout << MPIManager::thread() << " : "
                    << addLocalAssemble.locked << " " << addLocalAssemble.notLocked << " "
                    << addLocalAssemble.timesLocked << std::endl;
          #endif
        }

        template< class GridFunction, class JacobianOperator, class Iterators, class IntegrandsTuple>
        void assemble ( const GridFunction &u, JacobianOperator &jOp, const Iterators& iterators, const IntegrandsTuple& integrandsTuple ) const
        {
          AddLocalAssemble<JacobianOperator> addLocalAssemble(jOp);
          assemble( u, jOp, iterators, integrandsTuple, addLocalAssemble, 10 );
        }

        // accessors
        const GridPartType &gridPart () const { return gridPart_; }

        std::size_t gridSizeInterior () const { return gridSizeInterior_; }

      protected:
        const GridPartType &gridPart_;
        mutable std::size_t gridSizeInterior_;
      };


      template <class GalerkinOperator >
      static std::size_t accumulateGridSize( const ThreadSafeValue< GalerkinOperator >& ops )
      {
        std::size_t s = ops.size();
        std::size_t sum = 0;
        for( std::size_t i=0; i<s; ++i )
          sum += ops[ i ].gridSizeInterior();
        return sum;
      }

    } // namespace Impl

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////


    // GalerkinOperator
    // ----------------

    template< class Integrands, class DomainFunction, class RangeFunction = DomainFunction >
    struct GalerkinOperator
      : public virtual Operator< DomainFunction, RangeFunction >
    {
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;

      typedef typename RangeFunctionType::GridPartType GridPartType;

      typedef Impl::LocalGalerkinOperator< Integrands >  LocalGalerkinOperatorImplType;
      typedef Impl::GalerkinOperator< GridPartType  >    GalerkinOperatorImplType;

      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      typedef ThreadIterator< GridPartType > ThreadIteratorType;

      template< class... Args >
      explicit GalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : iterators_( gridPart ),
          opImpl_( gridPart ),
          localOp_( gridPart, std::forward< Args >( args )... ),
          gridSizeInterior_( 0 ),
          communicate_( true )
      {
      }

      void setCommunicate( const bool communicate )
      {
        communicate_ = communicate;
        if( ! communicate_ && Dune::Fem::Parameter::verbose() )
        {
          std::cout << "GalerkinOperator::setCommunicate: communicate was disabled!" << std::endl;
        }
      }

      void setQuadratureOrders(unsigned int interior, unsigned int surface)
      {
        size_t size = localOp_.size();
        for( size_t i=0; i<size; ++i )
          localOp_[ i ].setQuadratureOrders(interior,surface);
      }

      virtual bool nonlinear() const final override
      {
        return localOperator().nonlinear();
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        evaluate( u, w );
      }

      const GridPartType &gridPart () const { return op().gridPart(); }

      typedef Integrands ModelType;
      typedef Integrands DirichletModelType;
      ModelType &model() const { return localOperator().model(); }

      [[deprecated("Use localOperator instead!")]]
      const LocalGalerkinOperatorImplType& impl() const { return localOperator(); }

      //! return local operator holding instance of integrands
      const LocalGalerkinOperatorImplType& localOperator() const { return *localOp_; }

      std::size_t gridSizeInterior () const { return gridSizeInterior_; }

    protected:
      //! for implementation purposes
      const GalerkinOperatorImplType& op() const { return *opImpl_; }

      template < class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        iterators_.update();
        w.clear();

        std::shared_mutex mutex;

        auto doEval = [this, &u, &w, &mutex] ()
        {
          // TODO: Move this to be a class variable
          std::tuple< const LocalGalerkinOperatorImplType& > integrands( localOperator() );
          this->op().evaluate( u, w, this->iterators_, integrands, mutex );
        };

        try {
          // execute in parallel
          MPIManager :: run ( doEval );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = Impl::accumulateGridSize( opImpl_ );
        }
        catch ( const SingleThreadModeError& e )
        {
          // reset w from previous entries
          w.clear();
          // re-run in single thread mode if previous attempt failed
          std::tuple< const LocalGalerkinOperatorImplType& > integrands( localOperator() );
          op().evaluate( u, w, iterators_, integrands );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = op().gridSizeInterior();
        }

        // synchronize result
        if( communicate_ )
          w.communicate();
      }

      mutable ThreadIteratorType iterators_;
      ThreadSafeValue< GalerkinOperatorImplType > opImpl_;
      ThreadSafeValue< LocalGalerkinOperatorImplType > localOp_;

      mutable std::size_t gridSizeInterior_;
      bool communicate_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Integrands, class JacobianOperator >
    class DifferentiableGalerkinOperator
      : public GalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef GalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType > BaseType;

      typedef typename BaseType :: LocalGalerkinOperatorImplType  LocalGalerkinOperatorImplType;
    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType    DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType     RangeDiscreteFunctionSpaceType;

      typedef DiagonalAndNeighborStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType > DiagonalAndNeighborStencilType;
      typedef DiagonalStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType >            DiagonalStencilType;

      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit DifferentiableGalerkinOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                const RangeDiscreteFunctionSpaceType &rSpace,
                                                Args &&... args )
        : BaseType( rSpace.gridPart(), std::forward< Args >( args )... ),
          dSpace_(dSpace), rSpace_(rSpace),
          domainSpaceSequence_(dSpace.sequence()),
          rangeSpaceSequence_(rSpace.sequence()),
          stencilDAN_(), stencilD_()
      {
          if( hasSkeleton() )
            stencilDAN_.reset( new DiagonalAndNeighborStencilType( dSpace_, rSpace_ ) );
          else
            stencilD_.reset( new DiagonalStencilType( dSpace_, rSpace_ ) );
      }

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        assemble( u, jOp );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return rSpace_;
      }

      using BaseType::localOperator;
      using BaseType::nonlinear;

    protected:
      using BaseType::op;

      bool hasSkeleton() const
      {
        std::tuple< const LocalGalerkinOperatorImplType& > integrands( localOperator() );
        return op().hasSkeleton( integrands );
      }

      void prepare( JacobianOperatorType& jOp ) const
      {
        if ( domainSpaceSequence_ != domainSpace().sequence()
             || rangeSpaceSequence_ != rangeSpace().sequence() )
        {
          domainSpaceSequence_ = domainSpace().sequence();
          rangeSpaceSequence_ = rangeSpace().sequence();
          if( hasSkeleton() )
          {
            assert( stencilDAN_ );
            stencilDAN_->update();
          }
          else
          {
            assert( stencilD_ );
            stencilD_->update();
          }
        }
        if( hasSkeleton() )
          jOp.reserve( *stencilDAN_ );
        else
          jOp.reserve( *stencilD_ );
        // set all entries to zero
        jOp.clear();
      }

      template < class GridFunction >
      void assemble( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        // reserve memory and clear entries
        {
          prepare( jOp );
          iterators_.update();
        }

        std::shared_mutex mutex;

        auto doAssemble = [this, &u, &jOp, &mutex] ()
        {
          std::tuple< const LocalGalerkinOperatorImplType& > integrands( localOperator() );
          this->op().assemble( u, jOp, this->iterators_, integrands, mutex );
        };

        try {
          // execute in parallel
          MPIManager :: run ( doAssemble );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = Impl::accumulateGridSize( this->opImpl_ );
        }
        catch ( const SingleThreadModeError& e )
        {
          // redo assemble since it failed previously
          jOp.clear();
          std::tuple< const LocalGalerkinOperatorImplType& > integrands( localOperator() );
          op().assemble( u, jOp, iterators_, integrands );
          // update number of interior elements as sum over threads
          gridSizeInterior_ = op().gridSizeInterior();
        }

        // note: assembly done without local contributions so need
        // to call flush assembly
        jOp.flushAssembly();
      }

      using BaseType::iterators_;
      using BaseType::gridSizeInterior_;

      const DomainDiscreteFunctionSpaceType &dSpace_;
      const RangeDiscreteFunctionSpaceType &rSpace_;

      mutable int domainSpaceSequence_, rangeSpaceSequence_;

      mutable std::unique_ptr< DiagonalAndNeighborStencilType >  stencilDAN_;
      mutable std::unique_ptr< DiagonalStencilType > stencilD_;
    };



    // AutomaticDifferenceGalerkinOperator
    // -----------------------------------

    template< class Integrands, class DomainFunction, class RangeFunction >
    class AutomaticDifferenceGalerkinOperator
      : public GalerkinOperator< Integrands, DomainFunction, RangeFunction >,
        public AutomaticDifferenceOperator< DomainFunction, RangeFunction >
    {
      typedef GalerkinOperator< Integrands, DomainFunction, RangeFunction > BaseType;
      typedef AutomaticDifferenceOperator< DomainFunction, RangeFunction > AutomaticDifferenceOperatorType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit AutomaticDifferenceGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : BaseType( gridPart, std::forward< Args >( args )... ), AutomaticDifferenceOperatorType()
      {}
    };



    // ModelDifferentiableGalerkinOperator
    // -----------------------------------

    template < class LinearOperator, class ModelIntegrands >
    struct ModelDifferentiableGalerkinOperator
      : public DifferentiableGalerkinOperator< ModelIntegrands, LinearOperator >
    {
      typedef DifferentiableGalerkinOperator< ModelIntegrands, LinearOperator > BaseType;

      typedef typename ModelIntegrands::ModelType ModelType;

      typedef typename LinearOperator::DomainFunctionType RangeFunctionType;
      typedef typename LinearOperator::RangeSpaceType DiscreteFunctionSpaceType;

      ModelDifferentiableGalerkinOperator ( ModelType &model, const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( dfSpace.gridPart(), model )
      {}

      template< class GridFunction >
      void apply ( const GridFunction &u, RangeFunctionType &w ) const
      {
        (*this)( u, w );
      }

      template< class GridFunction >
      void apply ( const GridFunction &u, LinearOperator &jOp ) const
      {
        (*this).jacobian( u, jOp );
      }
    };

    namespace Impl
    {

      // GalerkinSchemeImpl
      // ------------------
      template< class Integrands, class LinearOperator, bool addDirichletBC,
                template <class,class> class DifferentiableGalerkinOperatorImpl >
      struct GalerkinSchemeTraits
      {
        template <class O, bool addDBC>
        struct DirichletBlockSelector { using type = void; };
        template <class O>
        struct DirichletBlockSelector<O,true> { using type = typename O::DirichletBlockVector; };

        using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
           DirichletWrapperOperator< DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
           DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator > >;
        using DirichletBlockVector = typename DirichletBlockSelector<
                 DirichletWrapperOperator<
                    DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
                 addDirichletBC>::type;

        typedef DifferentiableOperatorType type;
      };

      template< class Integrands, class LinearOperator, class LinearInverseOperator, bool addDirichletBC,
                template <class,class> class DifferentiableGalerkinOperatorImpl = DifferentiableGalerkinOperator >
      struct GalerkinSchemeImpl : public FemScheme< typename
                                  GalerkinSchemeTraits< Integrands, LinearOperator,
                                       addDirichletBC, DifferentiableGalerkinOperatorImpl>::type, // Operator
                                      LinearInverseOperator > // LinearInverseOperator
      {
        typedef FemScheme< typename GalerkinSchemeTraits< Integrands, LinearOperator,
                               addDirichletBC, DifferentiableGalerkinOperatorImpl>::type, // Operator
                               LinearInverseOperator > // LinearInverseOperator
                        BaseType;

        typedef typename BaseType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;

        GalerkinSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
                             const Integrands &integrands,
                             const ParameterReader& parameter = Parameter::container() )
          : BaseType(dfSpace,
                     parameter,
                     std::move(integrands))
        {}
      };

    } // end namespace Impl

    // GalerkinScheme
    // --------------

    template< class Integrands, class LinearOperator, class InverseOperator, bool addDirichletBC >
    using GalerkinScheme = Impl::GalerkinSchemeImpl< Integrands, LinearOperator, InverseOperator, addDirichletBC,
                                                     DifferentiableGalerkinOperator >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
