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
#include <dune/fem/schemes/integrands.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

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

      private:
        typedef CachingQuadrature< GridPartType, 0 > InteriorQuadratureType;
        typedef CachingQuadrature< GridPartType, 1 > SurfaceQuadratureType;

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
          InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder(w.order()) );

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

          InteriorQuadratureType quadrature( u.entity(), interiorQuadratureOrder( std::max(domainBasis.order(),rangeBasis.order()) ) );
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
          for( const auto qp : SurfaceQuadratureType( gridPart(), intersection, surfaceQuadratureOrder(w.order()), SurfaceQuadratureType::INSIDE ) )
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

          for( const auto qp : SurfaceQuadratureType( gridPart(), intersection, surfaceQuadratureOrder(rangeBasis.order()), SurfaceQuadratureType::INSIDE ) )
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
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, surfaceQuadratureOrder(wIn.order()), false );
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
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, surfaceQuadratureOrder(std::max( wIn.order(), wOut.order() )), false );
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
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &jInIn, J &jOutIn ) const
        {
          const auto &domainBasisIn = jInIn.domainBasisFunctionSet();
          const auto &domainBasisOut = jOutIn.domainBasisFunctionSet();

          const auto &rangeBasisIn = jInIn.rangeBasisFunctionSet();

          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, surfaceQuadratureOrder(rangeBasisIn.order()), false );
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
        void addLinearizedSkeletonIntegral ( const Intersection &intersection, const U &uIn, const U &uOut, DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &jInIn, J &jOutIn, J &jInOut, J &jOutOut ) const
        {
          const auto &domainBasisIn = jInIn.domainBasisFunctionSet();
          const auto &domainBasisOut = jOutIn.domainBasisFunctionSet();

          const auto &rangeBasisIn = jInIn.rangeBasisFunctionSet();
          const auto &rangeBasisOut = jInOut.rangeBasisFunctionSet();

          const auto geometry = intersection.geometry();
          const IntersectionQuadrature< SurfaceQuadratureType, conforming > quadrature( gridPart(), intersection, surfaceQuadratureOrder(std::max( rangeBasisIn.order(), rangeBasisOut.order() )), false );
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
          : gridPart_( gridPart ), integrands_( std::forward< Args >( args )... ),
            interiorQuadOrder_(0), surfaceQuadOrder_(0)
        {}
        void setQuadratureOrders(unsigned int interior, unsigned int surface)
        {
          interiorQuadOrder_ = interior;
          surfaceQuadOrder_ = surface;
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

          TemporaryLocalFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > wLocal( w.space() );

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

          TemporaryLocalFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > wInside( w.space() ), wOutside( w.space() );

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
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename DiscreteFunction::GridPartType, GridPartType >::value, "Argument 'w' and Integrands must be defined on the same grid part." );

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
          typedef TemporaryLocalMatrix< typename JacobianOperator::DomainSpaceType, typename JacobianOperator::RangeSpaceType > TemporaryLocalMatrixType;

          DiagonalStencil< typename JacobianOperator::DomainSpaceType, typename JacobianOperator::RangeSpaceType > stencil( jOp.domainSpace(), jOp.rangeSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;
          DomainValueVectorType phi = makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpLocal( jOp.domainSpace(), jOp.rangeSpace() );

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
          typedef TemporaryLocalMatrix< typename JacobianOperator::DomainSpaceType, typename JacobianOperator::RangeSpaceType > TemporaryLocalMatrixType;

          DiagonalAndNeighborStencil< typename JacobianOperator::DomainSpaceType, typename JacobianOperator::RangeSpaceType > stencil( jOp.domainSpace(), jOp.rangeSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;
          DomainValueVectorType phiIn = makeDomainValueVector( maxNumLocalDofs );
          DomainValueVectorType phiOut = makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpInIn( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutIn( jOp.domainSpace(), jOp.rangeSpace() );
          TemporaryLocalMatrixType jOpInOut( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutOut( jOp.domainSpace(), jOp.rangeSpace() );

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
          static_assert( std::is_same< typename GridFunction::GridPartType, GridPartType >::value, "Argument 'u' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::DomainSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );
          static_assert( std::is_same< typename JacobianOperator::RangeSpaceType::GridPartType, GridPartType >::value, "Argument 'jOp' and Integrands must be defined on the same grid part." );

          if( integrands_.hasSkeleton() )
            assemble( u, jOp, std::true_type() );
          else
            assemble( u, jOp, std::false_type() );
        }

        // accessors

        const GridPartType &gridPart () const { return gridPart_; }

        unsigned int interiorQuadratureOrder(unsigned int order) const { return interiorQuadOrder_==0 ? 2*order+3:interiorQuadOrder_; }
        unsigned int surfaceQuadratureOrder(unsigned int order)  const { return surfaceQuadOrder_==0 ? 2*order+3:surfaceQuadOrder_; }

      private:
        const GridPartType &gridPart_;
        mutable IntegrandsType integrands_;
        unsigned int interiorQuadOrder_;
        unsigned int surfaceQuadOrder_;

        mutable std::vector< RangeValueVectorType > rangeValues_;
        mutable RangeValueVectorType  values_;
        mutable DomainValueVectorType basisValues_;
        mutable DomainValueVectorType domainValues_;

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
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit DifferentiableGalerkinOperator ( const DomainSpaceType &dSpace, const RangeSpaceType &rSpace,
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

      const DomainSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeSpaceType& rangeSpace() const
      {
        return rSpace_;
      }

    protected:
      using BaseType::impl_;
      const DomainSpaceType &dSpace_;
      const RangeSpaceType &rSpace_;
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



    // GalerkinScheme
    // --------------

    template< class Integrands, class LinearOperator, class InverseOperator, bool addDirichletBC >
    struct GalerkinScheme
    {
      typedef InverseOperator InverseOperatorType;
      typedef Integrands ModelType;
      using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
         DirichletWrapperOperator< DifferentiableGalerkinOperator< Integrands, LinearOperator >>,
         DifferentiableGalerkinOperator< Integrands, LinearOperator > >;

      typedef typename DifferentiableOperatorType::DomainFunctionType DomainFunctionType;
      typedef typename DifferentiableOperatorType::RangeFunctionType RangeFunctionType;
      typedef typename DifferentiableOperatorType::JacobianOperatorType LinearOperatorType;

      typedef RangeFunctionType DiscreteFunctionType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef LinearOperator JacobianOperatorType;

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

      GalerkinScheme ( const DiscreteFunctionSpaceType &dfSpace, Integrands &integrands, ParameterReader parameter = Parameter::container() )
        : dfSpace_( dfSpace ),
          fullOperator_( dfSpace, dfSpace, std::move( integrands ) ),
          parameter_( std::move( parameter ) ),
          linearOperator_( "assembled elliptic operator", dfSpace, dfSpace ),
          invOp_(parameter_)
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
        invOp_.bind(fullOperator());
        DiscreteFunctionType rhs0 = rhs;
        setZeroConstraints( rhs0 );
        invOp_( rhs0, solution );
        invOp_.unbind();
        return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations() );
      }
      SolverInfo solve ( DiscreteFunctionType &solution ) const
      {
        DiscreteFunctionType bnd( solution );
        bnd.clear();
        return solve(bnd, solution);
      }

      template< class GridFunction >
      const LinearOperatorType &assemble ( const GridFunction &ubar )
      {
        fullOperator().jacobian( ubar, linearOperator_ );
        return linearOperator_;
      }

      const DiscreteFunctionSpaceType &space () const { return dfSpace_; }
      const GridPartType &gridPart () const { return space().gridPart(); }
      ModelType &model() const { return fullOperator().model(); }

      std::enable_if_t<addDirichletBC,void>
      setConstraints( DomainFunctionType &u ) const
      {
        fullOperator().setConstraints( u );
      }
      std::enable_if_t<addDirichletBC,void>
      setConstraints( const typename DiscreteFunctionType::RangeType &value, DiscreteFunctionType &u ) const
      {
        fullOperator().setConstraints( value, u );
      }
      std::enable_if_t<addDirichletBC,void>
      setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        fullOperator().setConstraints( u, v );
      }
      std::enable_if_t<addDirichletBC,void>
      subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
      {
        fullOperator().subConstraints( u, v );
      }
    protected:
      std::enable_if_t<addDirichletBC,void>
      setZeroConstraints( DiscreteFunctionType &u ) const { fullOperator().setConstraints( typename DiscreteFunctionType::RangeType(0), u ); }
      void setZeroConstraints( ... ) const { }
      const DiscreteFunctionSpaceType &dfSpace_;
      DifferentiableOperatorType fullOperator_;
      ParameterReader parameter_;
      LinearOperatorType linearOperator_;
      mutable NewtonOperatorType invOp_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_GALERKIN_HH
