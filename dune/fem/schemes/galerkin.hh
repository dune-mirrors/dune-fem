#ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
#define DUNE_FEM_SCHEMES_GALERKIN_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/hybridutilities.hh>

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
#include <dune/fem/misc/l2norm.hh>

#include <dune/fem/operator/common/localmatrixcolumn.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/schemes/integrands.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>

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
        // check for 'int order () const' or 'int order()'
        // and variants returning unsigned int or size_t
        template <class T> static std::true_type testSignature(int (T::*)() const);
        template <class T> static std::true_type testSignature(int (T::*)());
        template <class T> static std::true_type testSignature(unsigned int (T::*)() const);
        template <class T> static std::true_type testSignature(unsigned int (T::*)());
        template <class T> static std::true_type testSignature(std::size_t (T::*)() const);
        template <class T> static std::true_type testSignature(std::size_t (T::*)());

        template <class T>
        static decltype(testSignature(&T::order)) test(std::nullptr_t);

        template <class T>
        static std::false_type test(...);

        using type = decltype(test<M>(nullptr));

        template <class F>
        static int callOrder(const F& f, std::false_type)
        {
#ifndef NDEBUG
          std::cerr << "WARNING: not order method available on " << typeid(F).name() << ", defaulting to 1!" << std::endl;
#endif
          return 1;
        }

        template <class F>
        static int callOrder(const F& f, std::true_type)
        {
          return f.order();
        }

      public:
        template <class F>
        static int order (const F& f ) { return callOrder(f, type() ); }
      };

      // GalerkinOperator
      // ----------------

      template< class Integrands >
      struct GalerkinOperator
      {
        typedef GalerkinOperator<Integrands> ThisType;
        typedef std::conditional_t< Fem::IntegrandsTraits< Integrands >::isFull, Integrands, FullIntegrands< Integrands > > IntegrandsType;

        typedef typename IntegrandsType::GridPartType GridPartType;

        typedef typename GridPartType::ctype ctype;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        template <class Space>
        struct QuadratureSelector
        {
          typedef CachingQuadrature< GridPartType, 0, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > InteriorQuadratureType;
          typedef CachingQuadrature< GridPartType, 1, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > SurfaceQuadratureType;
        // typedef CachingQuadrature< GridPartType, 0, Dune::FemPy::FempyQuadratureTraits > InteriorQuadratureType;
        // typedef CachingQuadrature< GridPartType, 1, Dune::FemPy::FempyQuadratureTraits > SurfaceQuadratureType;
        };

      private:
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
          Hybrid::forEach( std::index_sequence_for< T... >(), [ &u, &x, &phi ] ( auto i ) { GalerkinOperator::value( u, x, std::get< i >( phi ) ); } );
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
          Hybrid::forEach( std::index_sequence_for< T... >(), [ &basis, &x, &phi ] ( auto i ) { GalerkinOperator::values( basis, x, std::get< i >( phi ) ); } );
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
          if( !integrands_.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();

          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: InteriorQuadratureType  InteriorQuadratureType;
          InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder(maxOrder(u, w)) );

          // evaluate u for all quadrature points
          DomainValueVectorType& domains = domainValues_;
          domainValue( u, quadrature, domains );

          auto& ranges = values_;
          resizeRangeValueVector( ranges, quadrature.nop() );

          // evaluate integrands for all quadrature points
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.position() );
            assignRange( ranges, qp.index(), integrands_.interior( qp, domainValue( qp.index(), domains ) ), weight );
          }

          // add to w for all quadrature points
          axpyQuadrature( w, quadrature, ranges );
        }

        template< class U, class J >
        void addLinearizedInteriorIntegral ( const U &u, DomainValueVectorType &phi, J &j ) const
        {
          if( !integrands_.init( u.entity() ) )
            return;

          const auto geometry = u.entity().geometry();
          const auto &domainBasis = j.domainBasisFunctionSet();
          const auto &rangeBasis = j.rangeBasisFunctionSet();

          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: InteriorQuadratureType  InteriorQuadratureType;
          InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder( maxOrder( u, domainBasis, rangeBasis )) );
          const size_t domainSize = domainBasis.size();
          const size_t quadNop = quadrature.nop();

          auto& basisValues = basisValues_;
          resizeDomainValueVector( basisValues_, domainSize );

          // evaluate u for all quadrature points
          DomainValueVectorType& domains = domainValues_;
          domainValue( u, quadrature, domains );

          rangeValues_.resize( domainSize );
          for( std::size_t col = 0; col < domainSize; ++col )
          {
            resizeRangeValueVector( rangeValues_[ col ], quadNop );
          }

          // evaluate all basis functions and integrands
          for( const auto qp : quadrature )
          {
            values( domainBasis, qp, basisValues );
            const auto weight = qp.weight() * geometry.integrationElement( qp.position() );
            auto integrand = integrands_.linearizedInterior( qp, domainValue( qp.index(), domains ) );
            for( std::size_t col = 0; col < domainSize; ++col )
            {
              assignRange( rangeValues_[ col ], qp.index(), integrand( value( basisValues, col ) ), weight );
            }
          }

          // add to local matrix for all quadrature points and basis functions
          for( std::size_t col = 0; col < domainSize; ++col )
          {
            LocalMatrixColumn< J > jCol( j, col );
            axpyQuadrature( jCol, quadrature, rangeValues_[ col ] );
          }
        }

        // boundary integral

        template< class Intersection, class U, class W >
        void addBoundaryIntegral ( const Intersection &intersection, const U &u, W &w ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          typedef typename QuadratureSelector< typename W::DiscreteFunctionSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          SurfaceQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder( u, w )), SurfaceQuadratureType::INSIDE );
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            RangeValueType integrand = integrands_.boundary( qp, domainValue( u, qp ) );

            Hybrid::forEach( RangeValueIndices(), [ &qp, &w, &integrand, weight ] ( auto i ) {
                std::get< i >( integrand ) *= weight;
                w.axpy( qp, std::get< i >( integrand ) );
              } );
          }
        }

        template< class Intersection, class U, class J >
        void addLinearizedBoundaryIntegral ( const Intersection &intersection, const U &u, DomainValueVectorType &phi, J &j ) const
        {
          if( !integrands_.init( intersection ) )
            return;

          const auto geometry = intersection.geometry();
          const auto &domainBasis = j.domainBasisFunctionSet();
          const auto &rangeBasis = j.rangeBasisFunctionSet();

          typedef typename QuadratureSelector< typename J::RangeSpaceType > :: SurfaceQuadratureType  SurfaceQuadratureType;
          SurfaceQuadratureType quadrature( gridPart(), intersection, surfaceQuadratureOrder(maxOrder(u, domainBasis, rangeBasis )), SurfaceQuadratureType::INSIDE );
          for( const auto qp : quadrature )
          {
            const ctype weight = qp.weight() * geometry.integrationElement( qp.localPosition() );

            values( domainBasis, qp, phi );
            auto integrand = integrands_.linearizedBoundary( qp, domainValue( u, qp ) );

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
        }

        // addSkeletonIntegral

      private:
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
            std::pair< RangeValueType, RangeValueType > integrand = integrands_.skeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );

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
            std::pair< RangeValueType, RangeValueType > integrand = integrands_.skeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );

            Hybrid::forEach( RangeValueIndices(), [ &qpIn, &wIn, &qpOut, &wOut, &integrand, weight ] ( auto i ) {
                std::get< i >( integrand.first ) *= weight;
                wIn.axpy( qpIn, std::get< i >( integrand.first ) );

                std::get< i >( integrand.second ) *= weight;
                wOut.axpy( qpOut, std::get< i >( integrand.second ) );
              } );
          }
        }

        template< bool conforming, class Intersection, class U, class J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut,
                                             DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &jInIn, J &jOutIn ) const
        {
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

            auto integrand = integrands_.linearizedSkeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );
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
                                             DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &jInIn, J &jOutIn, J &jInOut, J &jOutOut ) const
        {
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

            auto integrand = integrands_.linearizedSkeleton( qpIn, domainValue( uIn, qpIn ), qpOut, domainValue( uOut, qpOut ) );
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
          if( !integrands_.init( intersection ) )
            return;

          if( intersection.conforming() )
            addSkeletonIntegral< true >( intersection, uIn, uOut, w... );
          else
            addSkeletonIntegral< false >( intersection, uIn, uOut, w... );
        }

        template< class Intersection, class U, class... J >
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &... j ) const
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
        explicit GalerkinOperator ( const GridPartType &gridPart, Args &&... args )
          : gridPart_( gridPart ),
            integrands_( std::forward< Args >( args )... ),
            defaultInteriorOrder_( [] (const int order) { return 2 * order; } ),
            defaultSurfaceOrder_ ( [] (const int order) { return 2 * order + 1; } ),
            interiorQuadOrder_(0), surfaceQuadOrder_(0),
            communicate_( true )
        {
        }

        void setCommunicate( const bool communicate )
        {
          communicate_ = communicate;
          if( ! communicate_ && Dune::Fem::Parameter::verbose() )
          {
            std::cout << "Impl::GalerkinOperator::setCommunicate: communicate was disabled!" << std::endl;
          }
        }

        void setQuadratureOrders(unsigned int interior, unsigned int surface)
        {
          interiorQuadOrder_ = interior;
          surfaceQuadOrder_  = surface;
        }

        IntegrandsType &model() const
        {
          return integrands_;
        }

      private:
        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, std::false_type ) const
        {
          w.clear();

          typedef typename DiscreteFunction::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wLocal( w.space() );

          defaultInteriorOrder_ = [] (const int order) { return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::volumeOrder(order); };
          defaultSurfaceOrder_  = [] (const int order) { return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::surfaceOrder(order); };

          Dune::Fem::ConstLocalFunction< GridFunction > uLocal( u );
          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            uLocal.bind( entity );
            wLocal.bind( entity );
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
        }

        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w, std::true_type ) const
        {
          w.clear();

          typedef typename DiscreteFunction::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

          TemporaryLocalFunction< DiscreteFunctionSpaceType > wInside( w.space() ), wOutside( w.space() );

          defaultInteriorOrder_ = [] (const int order) { return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::volumeOrder(order); };
          defaultSurfaceOrder_  = [] (const int order) { return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::surfaceOrder(order); };

          Dune::Fem::ConstLocalFunction< GridFunction > uInside( u );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            uInside.bind( inside );
            wInside.bind( inside );
            wInside.clear();

            if( integrands_.hasInterior() )
              addInteriorIntegral( uInside, wInside );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              // check neighbor first since on periodic boundaries both,
              // neighbor and boundary are true, so we treat neighbor first
              if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                Dune::Fem::ConstLocalFunction< GridFunction > uOutside( u );
                uOutside.bind( outside );
                if( outside.partitionType() != InteriorEntity )
                  addSkeletonIntegral( intersection, uInside, uOutside, wInside );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  wOutside.bind( outside );
                  wOutside.clear();
                  addSkeletonIntegral( intersection, uInside, uOutside, wInside, wOutside );
                  w.addLocalDofs( outside, wOutside.localDofVector() );
                }
              }
              else if( intersection.boundary() )
              {
                if( integrands_.hasBoundary() )
                  addBoundaryIntegral( intersection, uInside, wInside );
              }
            }

            w.addLocalDofs( inside, wInside.localDofVector() );
          }
        }

      public:
        template< class GridFunction, class DiscreteFunction >
        void evaluate ( const GridFunction &u, DiscreteFunction &w ) const
        {
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename DiscreteFunction::GridPartType, GridPartType >::value, "Argument 'w' and Integrands must be defined on the same grid part." );

          if( integrands_.hasSkeleton() )
            evaluate( u, w, std::true_type() );
          else
            evaluate( u, w, std::false_type() );

          // synchronize result
          if( communicate_ )
            w.communicate();
        }

        // assemble

      private:
        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp, std::false_type ) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;

          // select correct default quadrature orders
          defaultInteriorOrder_ = [] (const int order) { return Capabilities::DefaultQuadrature< RangeSpaceType >::volumeOrder(order); };
          defaultSurfaceOrder_  = [] (const int order) { return Capabilities::DefaultQuadrature< RangeSpaceType >::surfaceOrder(order); };

          DiagonalStencil< DomainSpaceType, RangeSpaceType > stencil( jOp.domainSpace(), jOp.rangeSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;
          DomainValueVectorType phi = makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpLocal( jOp.domainSpace(), jOp.rangeSpace() );
          Dune::Fem::ConstLocalFunction< GridFunction > uLocal( u );

          for( const EntityType &entity : elements( gridPart(), Partitions::interiorBorder ) )
          {
            uLocal.bind( entity );

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
        }

        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp, std::true_type ) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;

          // select correct default quadrature orders
          defaultInteriorOrder_ = [] (const int order) { return Capabilities::DefaultQuadrature< RangeSpaceType >::volumeOrder(order); };
          defaultSurfaceOrder_  = [] (const int order) { return Capabilities::DefaultQuadrature< RangeSpaceType >::surfaceOrder(order); };

          DiagonalAndNeighborStencil< DomainSpaceType, RangeSpaceType > stencil( jOp.domainSpace(), jOp.rangeSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;
          DomainValueVectorType phiIn = makeDomainValueVector( maxNumLocalDofs );
          DomainValueVectorType phiOut = makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpInIn( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutIn( jOp.domainSpace(), jOp.rangeSpace() );
          TemporaryLocalMatrixType jOpInOut( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutOut( jOp.domainSpace(), jOp.rangeSpace() );
          Dune::Fem::ConstLocalFunction< GridFunction > uIn( u );

          const auto &indexSet = gridPart().indexSet();
          for( const EntityType &inside : elements( gridPart(), Partitions::interiorBorder ) )
          {
            uIn.bind( inside );

            jOpInIn.init( inside, inside );
            jOpInIn.clear();

            if( integrands_.hasInterior() )
              addLinearizedInteriorIntegral( uIn, phiIn, jOpInIn );

            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              // check neighbor first since on periodic boundaries both,
              // neighbor and boundary are true, so we treat neighbor first
              if( intersection.neighbor() )
              {
                const EntityType &outside = intersection.outside();

                jOpOutIn.init( outside, inside );
                jOpOutIn.clear();

                Dune::Fem::ConstLocalFunction< GridFunction > uOut( u );
                uOut.bind( outside );

                if( outside.partitionType() != InteriorEntity )
                  addLinearizedSkeletonIntegral( intersection, uIn, uOut, phiIn, phiOut, jOpInIn, jOpOutIn );
                else if( indexSet.index( inside ) < indexSet.index( outside ) )
                {
                  jOpInOut.init( inside, outside );
                  jOpInOut.clear();
                  jOpOutOut.init( outside, outside );
                  jOpOutOut.clear();

                  addLinearizedSkeletonIntegral( intersection, uIn, uOut, phiIn, phiOut, jOpInIn, jOpOutIn, jOpInOut, jOpOutOut );

                  jOp.addLocalMatrix( inside, outside, jOpInOut );
                  jOp.addLocalMatrix( outside, outside, jOpOutOut );
                }

                jOp.addLocalMatrix( outside, inside, jOpOutIn );
              }
              else if( intersection.boundary() )
              {
                if( integrands_.hasBoundary() )
                  addLinearizedBoundaryIntegral( intersection, uIn, phiIn, jOpInIn );
              }
            }

            jOp.addLocalMatrix( inside, inside, jOpInIn );
          }
        }

      public:
        template< class GridFunction, class JacobianOperator >
        void assemble ( const GridFunction &u, JacobianOperator &jOp ) const
        {
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::DomainSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::RangeSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );

          if( integrands_.hasSkeleton() )
            assemble( u, jOp, std::true_type() );
          else
            assemble( u, jOp, std::false_type() );

          // note: assembly done without local contributions so need
          // to call flush assembly
          jOp.flushAssembly();
        }

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
        mutable RangeValueVectorType  values_;
        mutable DomainValueVectorType basisValues_;
        mutable DomainValueVectorType domainValues_;

        bool communicate_;
      };

    } // namespace Impl




    // GalerkinOperator
    // ----------------

    template< class Integrands, class DomainFunction, class RangeFunction = DomainFunction >
    struct GalerkinOperator
      : public virtual Operator< DomainFunction, RangeFunction >
    {
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;

      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      typedef typename RangeFunctionType::GridPartType GridPartType;

      template< class... Args >
      explicit GalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : impl_( gridPart, std::forward< Args >( args )... )
      {}

      void setCommunicate( const bool communicate ) { impl_.setCommunicate( communicate ); }

      void setQuadratureOrders(unsigned int interior, unsigned int surface) { impl_.setQuadratureOrders(interior,surface); }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        impl_.evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        return impl_.evaluate( u, w );
      }

      const GridPartType &gridPart () const { return impl_.gridPart(); }

      typedef Integrands ModelType;
      typedef Integrands DirichletModelType;
      ModelType &model() const { return impl_.model(); }

    protected:
      Impl::GalerkinOperator< Integrands > impl_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Integrands, class JacobianOperator >
    class DifferentiableGalerkinOperator
      : public GalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType >,
        public DifferentiableOperator< JacobianOperator >
    {
      typedef GalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType > BaseType;

    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit DifferentiableGalerkinOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                const RangeDiscreteFunctionSpaceType &rSpace,
                                                Args &&... args )
        : BaseType( rSpace.gridPart(), std::forward< Args >( args )... ),
          dSpace_(dSpace), rSpace_(rSpace)
      {}

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        impl_.assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        impl_.assemble( u, jOp );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return rSpace_;
      }

    protected:
      using BaseType::impl_;
      const DomainDiscreteFunctionSpaceType &dSpace_;
      const RangeDiscreteFunctionSpaceType &rSpace_;
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
      template <class O, bool addDirichletBC>
      struct DirichletBlockSelector { using type = void; };
      template <class O>
      struct DirichletBlockSelector<O,true> { using type = typename O::DirichletBlockVector; };
      template< class Integrands, class LinearOperator, class InverseOperator, bool addDirichletBC,
                template <class,class> class DifferentiableGalerkinOperatorImpl = DifferentiableGalerkinOperator >
      struct GalerkinSchemeImpl
      {
        typedef InverseOperator InverseOperatorType;
        typedef Integrands ModelType;
        using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
           DirichletWrapperOperator< DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
           DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator > >;
        using DirichletBlockVector = typename DirichletBlockSelector<
                 DirichletWrapperOperator<
                    DifferentiableGalerkinOperatorImpl< Integrands, LinearOperator >>,
                 addDirichletBC>::type;

        typedef typename DifferentiableOperatorType::DomainFunctionType DomainFunctionType;
        typedef typename DifferentiableOperatorType::RangeFunctionType RangeFunctionType;
        typedef typename DifferentiableOperatorType::JacobianOperatorType LinearOperatorType;
        typedef typename DifferentiableOperatorType::JacobianOperatorType JacobianOperatorType;

        typedef RangeFunctionType DiscreteFunctionType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

        typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, InverseOperator > NewtonOperatorType;
        typedef InverseOperator LinearInverseOperatorType;
        typedef typename NewtonOperatorType::ErrorMeasureType ErrorMeasureType;

        struct SolverInfo
        {
          SolverInfo ( bool converged, int linearIterations, int nonlinearIterations )
            : converged( converged ), linearIterations( linearIterations ), nonlinearIterations( nonlinearIterations )
          {}

          bool converged;
          int linearIterations, nonlinearIterations;
        };

        GalerkinSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
                             const Integrands &integrands,
                             const ParameterReader& parameter = Parameter::container() )
          : dfSpace_( dfSpace ),
            fullOperator_( dfSpace, dfSpace, std::move( integrands ) ),
            invOp_(parameter)
        {}

        void setQuadratureOrders(unsigned int interior, unsigned int surface) { fullOperator().setQuadratureOrders(interior,surface); }

        const DifferentiableOperatorType &fullOperator() const { return fullOperator_; }
        DifferentiableOperatorType &fullOperator() { return fullOperator_; }

        void constraint ( DiscreteFunctionType &u ) const {}

        template< class GridFunction >
        void operator() ( const GridFunction &u, DiscreteFunctionType &w ) const
        {
          fullOperator()( u, w );
        }

        void setErrorMeasure(ErrorMeasureType &errorMeasure) const
        {
          invOp_.setErrorMeasure(errorMeasure);
        }

        SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution ) const
        {
          DiscreteFunctionType rhs0 = rhs;
          setZeroConstraints( rhs0 );
          setModelConstraints( solution );

          invOp_.bind(fullOperator());
          invOp_( rhs0, solution );
          invOp_.unbind();
          return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations() );
        }

        SolverInfo solve ( DiscreteFunctionType &solution ) const
        {
          DiscreteFunctionType bnd( solution );
          bnd.clear();
          setModelConstraints( solution );
          invOp_.bind(fullOperator());
          invOp_( bnd, solution );
          invOp_.unbind();
          return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations() );
        }

        template< class GridFunction >
        void jacobian( const GridFunction &ubar, LinearOperatorType &linearOp) const
        {
          fullOperator().jacobian( ubar, linearOp );
        }

        const DiscreteFunctionSpaceType &space () const { return dfSpace_; }
        const GridPartType &gridPart () const { return space().gridPart(); }
        ModelType &model() const { return fullOperator().model(); }

        void setConstraints( DomainFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u );
        }
        void setConstraints( const typename DiscreteFunctionType::RangeType &value, DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( value, u );
        }
        void setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u, v );
        }
        void subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().subConstraints( u, v );
        }
        const auto& dirichletBlocks() const
        {
          if constexpr (addDirichletBC)
            return fullOperator().dirichletBlocks();
        }

      protected:
        void setZeroConstraints( DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( typename DiscreteFunctionType::RangeType(0), u );
        }
        void setModelConstraints( DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u );
        }
        const DiscreteFunctionSpaceType &dfSpace_;
        DifferentiableOperatorType fullOperator_;
        mutable NewtonOperatorType invOp_;
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
