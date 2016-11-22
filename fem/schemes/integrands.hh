#ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
#define DUNE_FEM_SCHEMES_INTEGRANDS_HH

#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/fempy/quadrature/cachingpoint.hh>
#include <dune/fempy/quadrature/elementpoint.hh>

namespace Dune
{

  namespace Fem
  {

    // IntegrandsTraits
    // ----------------

    namespace Impl
    {

      namespace IntegrandsTraits
      {

        template< class Integrands >
        using ValueType = typename Integrands::ValueType;

        template< class Integrands >
        using GridPartType = typename Integrands::GridPartType;


        template< class Integrands >
        using EntityType = typename GridPartType< Integrands >::template Codim< 0 >::EntityType;

        template< class Integrands >
        using IntersectionType = typename GridPartType< Integrands >::IntersectionType;


        template< class Integrands >
        using InteriorQuadratureType = CachingQuadrature< GridPartType< Integrands >, 0 >;

        template< class Integrands >
        using SurfaceQuadratureType = CachingQuadrature< GridPartType< Integrands >, 1 >;


        template< class Integrands >
        using InteriorQuadraturePointType = QuadraturePointWrapper< InteriorQuadratureType< Integrands > >;

        template< class Integrands >
        using SurfaceQuadraturePointType = QuadraturePointWrapper< SurfaceQuadratureType< Integrands > >;


        template< class Integrands >
        static std::true_type interior ( const Integrands &, decltype( std::declval< const Integrands & >().interior( std::declval< const InteriorQuadraturePointType< Integrands > & >(), std::declval< const ValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type interior ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasInterior() ), bool >::value, int > = 0 >
        static std::true_type hasInterior ( const Integrands & );

        static std::false_type hasInterior ( ... );


        template< class Integrands >
        static std::true_type boundary ( const Integrands &, decltype( std::declval< const Integrands & >().boundary( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const ValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type boundary ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasBoundary() ), bool >::value, int > = 0 >
        static std::true_type hasBoundary ( const Integrands & );

        static std::false_type hasBoundary ( ... );


        template< class Integrands >
        static std::true_type skeleton ( const Integrands &, decltype( std::declval< const Integrands & >().skeleton( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const ValueType< Integrands > & >(), std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const ValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type skeleton ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasSkeleton() ), bool >::value, int > = 0 >
        static std::true_type hasSkeleton ( const Integrands & );

        static std::false_type hasSkeleton ( ... );

      } // namespace IntegrandsTraits

    } // namespace Impl

    template< class Integrands >
    struct IntegrandsTraits
    {
      typedef Impl::IntegrandsTraits::ValueType< Integrands > ValueType;
      typedef Impl::IntegrandsTraits::GridPartType< Integrands > GridPartType;

      typedef Impl::IntegrandsTraits::EntityType< Integrands > EntityType;
      typedef Impl::IntegrandsTraits::IntersectionType< Integrands > IntersectionType;

      static const bool interior = decltype( Impl::IntegrandsTraits::interior( std::declval< const Integrands & >() ) )::value;
      static const bool hasInterior = decltype( Impl::IntegrandsTraits::hasInterior( std::declval< const Integrands & >() ) )::value;

      static const bool boundary = decltype( Impl::IntegrandsTraits::boundary( std::declval< const Integrands & >() ) )::value;
      static const bool hasBoundary = decltype( Impl::IntegrandsTraits::hasBoundary( std::declval< const Integrands & >() ) )::value;

      static const bool skeleton = decltype( Impl::IntegrandsTraits::skeleton( std::declval< const Integrands & >() ) )::value;
      static const bool hasSkeleton = decltype( Impl::IntegrandsTraits::hasSkeleton( std::declval< const Integrands & >() ) )::value;

      static_assert( (!hasInterior || interior), "Existence of method 'hasInterior' implies existence of method interior." );
      static_assert( (!hasBoundary || boundary), "Existence of method 'hasBoundary' implies existence of method boundary." );
      static_assert( (!hasSkeleton || skeleton), "Existence of method 'hasSkeleton' implies existence of method skeleton." );

      static const bool isFull = hasInterior && hasBoundary && hasSkeleton;
    };



    // FullIntegrands
    // --------------

    template< class Integrands >
    struct FullIntegrands
    {
      typedef typename IntegrandsTraits< Integrands >::ValueType ValueType;
      typedef typename IntegrandsTraits< Integrands >::GridPartType GridPartType;

      typedef typename IntegrandsTraits< Integrands >::EntityType EntityType;
      typedef typename IntegrandsTraits< Integrands >::IntersectionType IntersectionType;

    private:
      template< class T, std::enable_if_t< IntegrandsTraits< T >::hasInterior, int > = 0 >
      static bool hasInterior ( const T &integrands )
      {
        return integrands.hasInterior();
      }

      template< class T, std::enable_if_t< !IntegrandsTraits< T >::hasInterior, int > = 0 >
      static bool hasInterior ( const T &integrands )
      {
        return IntegrandsTraits< T >::interior;
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::interior, int > = 0 >
      static ValueType interior ( const T &integrands, const Point &x, const ValueType &u )
      {
        return integrands.interior( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::interior, int > = 0 >
      static ValueType interior ( const T &integrands, const Point &x, const ValueType &u )
      {
        return ValueType();
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::interior, int > = 0 >
      static auto linearizedInterior ( const T &integrands, const Point &x, const ValueType &u )
      {
        return integrands.linearizedInterior( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::interior, int > = 0 >
      static auto linearizedInterior ( const T &integrands, const Point &x, const ValueType &u )
      {
        return [] ( const ValueType & ) { return ValueType(); };
      }

      template< class T, std::enable_if_t< IntegrandsTraits< T >::hasBoundary, int > = 0 >
      static bool hasBoundary ( const T &integrands )
      {
        return integrands.hasBoundary();
      }

      template< class T, std::enable_if_t< !IntegrandsTraits< T >::hasBoundary, int > = 0 >
      static bool hasBoundary ( const T &integrands )
      {
        return IntegrandsTraits< T >::boundary;
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::boundary, int > = 0 >
      static ValueType boundary ( const T &integrands, const Point &x, const ValueType &u )
      {
        return integrands.boundary( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::boundary, int > = 0 >
      static ValueType boundary ( const T &integrands, const Point &x, const ValueType &u )
      {
        return ValueType();
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::boundary, int > = 0 >
      static auto linearizedBoundary ( const T &integrands, const Point &x, const ValueType &u )
      {
        return integrands.linearizedBoundary( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::boundary, int > = 0 >
      static auto linearizedBoundary ( const T &integrands, const Point &x, const ValueType &u )
      {
        return [] ( const ValueType & ) { return ValueType(); };
      }

      template< class T, std::enable_if_t< IntegrandsTraits< T >::hasSkeleton, int > = 0 >
      static bool hasSkeleton ( const T &integrands )
      {
        return integrands.hasSkeleton();
      }

      template< class T, std::enable_if_t< !IntegrandsTraits< T >::hasSkeleton, int > = 0 >
      static bool hasSkeleton ( const T &integrands )
      {
        return IntegrandsTraits< T >::skeleton;
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::skeleton, int > = 0 >
      static std::pair< ValueType, ValueType > skeleton ( const T &integrands, const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut )
      {
        return integrands.skeleton( xIn, uIn, xOut, uOut );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::skeleton, int > = 0 >
      static std::pair< ValueType, ValueType > skeleton ( const T &integrands, const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut )
      {
        return std::make_pair( ValueType(), ValueType() );
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::skeleton, int > = 0 >
      static auto linearizedSkeleton ( const T &integrands, const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut )
      {
        return integrands.linearizedSkeleton( xIn, uIn, xOut, uOut );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::skeleton, int > = 0 >
      static auto linearizedSkeleton ( const T &integrands, const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut )
      {
        auto zero = [] ( const ValueType & ) { return std::make_pair( ValueType(), ValueType() ); };
        return std::make_pair( zero, zero );
      }

    public:
      template< class... Args >
      explicit FullIntegrands ( Args &&... args )
        : integrands_( std::forward< Args >( args )... )
      {}

      bool init ( const EntityType &entity ) { return std::ref( integrands_ ).get().init( entity ); }
      bool init ( const IntersectionType &intersection ) { return std::ref( integrands_ ).get().init( intersection ); }

      bool hasInterior () const { return hasInterior( std::ref( integrands_ ).get() ); }

      template< class Point >
      ValueType interior ( const Point &x, const ValueType &u ) const
      {
        return interior( std::ref( integrands_ ).get(), x, u );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return linearizedInterior( std::ref( integrands_ ).get(), x, u );
      }

      bool hasBoundary () const { return hasBoundary( std::ref( integrands_ ).get() ); }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        return boundary( std::ref( integrands_ ).get(), x, u );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return linearizedBoundary( std::ref( integrands_ ).get(), x, u );
      }

      bool hasSkeleton () const { return hasSkeleton( std::ref( integrands_ ).get() ); }

      template< class Point >
      std::pair< ValueType, ValueType > skeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        return skeleton( std::ref( integrands_ ).get(), xIn, uIn, xOut, uOut );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        return linearizedSkeleton( std::ref( integrands_ ).get(), xIn, uIn, xOut, uOut );
      }

    private:
      Integrands integrands_;
    };



    // VirtualizedIntegrands
    // ---------------------

    template< class GridPart, class Value >
    class VirtualizedIntegrands
    {
      typedef VirtualizedIntegrands< GridPart, Value > This;

    public:
      typedef GridPart GridPartType;
      typedef Value ValueType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

    private:
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      typedef FemPy::CachingPoint< LocalCoordinateType, 0 > InteriorCachingPointType;
      typedef FemPy::ElementPoint< LocalCoordinateType, 0 > InteriorElementPointType;
      typedef FemPy::CachingPoint< LocalCoordinateType, 1 > SurfaceCachingPointType;
      typedef FemPy::ElementPoint< LocalCoordinateType, 1 > SurfaceElementPointType;

      template< class QP >
      static Fem::QuadraturePointWrapper< QP > asQP ( const QP &qp )
      {
        return static_cast< Fem::QuadraturePointWrapper< QP > >( qp );
      }

      template< class R >
      using Linearization = std::function< R( const ValueType & ) >;

      template< class R >
      using LinearizationPair = std::pair< Linearization< std::pair< R, R > >, Linearization< std::pair< R, R > > >;

      struct Interface
      {
        virtual ~Interface ()  {}
        virtual Interface *clone () const = 0;

        virtual bool init ( const EntityType &entity ) = 0;
        virtual bool init ( const IntersectionType &intersection ) = 0;

        virtual bool hasInterior () const = 0;
        virtual ValueType interior ( const InteriorCachingPointType &x, const ValueType &u ) const = 0;
        virtual ValueType interior ( const InteriorElementPointType &x, const ValueType &u ) const = 0;
        virtual Linearization< ValueType > linearizedInterior ( const InteriorCachingPointType &x, const ValueType &u ) const = 0;
        virtual Linearization< ValueType > linearizedInterior ( const InteriorElementPointType &x, const ValueType &u ) const = 0;

        virtual bool hasBoundary () const = 0;
        virtual ValueType boundary ( const SurfaceCachingPointType &x, const ValueType &u ) const = 0;
        virtual ValueType boundary ( const SurfaceElementPointType &x, const ValueType &u ) const = 0;
        virtual Linearization< ValueType > linearizedBoundary ( const SurfaceCachingPointType &x, const ValueType &u ) const = 0;
        virtual Linearization< ValueType > linearizedBoundary ( const SurfaceElementPointType &x, const ValueType &u ) const = 0;

        virtual bool hasSkeleton () const = 0;
        virtual std::pair< ValueType, ValueType > skeleton ( const SurfaceCachingPointType &xIn, const ValueType &uIn, const SurfaceCachingPointType &xOut, const ValueType &uOut ) const = 0;
        virtual std::pair< ValueType, ValueType > skeleton ( const SurfaceElementPointType &xIn, const ValueType &uIn, const SurfaceElementPointType &xOut, const ValueType &uOut ) const = 0;
        virtual LinearizationPair< ValueType > linearizedSkeleton ( const SurfaceCachingPointType &xIn, const ValueType &uIn, const SurfaceCachingPointType &xOut, const ValueType &uOut ) const = 0;
        virtual LinearizationPair< ValueType > linearizedSkeleton ( const SurfaceElementPointType &xIn, const ValueType &uIn, const SurfaceElementPointType &xOut, const ValueType &uOut ) const = 0;
      };

      template< class Impl >
      struct Implementation final
        : public Interface
      {
        Implementation ( Impl impl ) : impl_( std::move( impl ) ) {}
        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual bool init ( const EntityType &entity ) override { return impl().init( entity ); }
        virtual bool init ( const IntersectionType &intersection ) override { return impl().init( intersection ); }

        virtual bool hasInterior () const override { return impl().hasInterior(); }
        virtual ValueType interior ( const InteriorCachingPointType &x, const ValueType &u ) const override { return impl().interior( asQP( x ), u ); }
        virtual ValueType interior ( const InteriorElementPointType &x, const ValueType &u ) const override { return impl().interior( asQP( x ), u ); }
        virtual Linearization< ValueType > linearizedInterior ( const InteriorCachingPointType &x, const ValueType &u ) const override { return impl().linearizedInterior( asQP( x ), u ); }
        virtual Linearization< ValueType > linearizedInterior ( const InteriorElementPointType &x, const ValueType &u ) const override { return impl().linearizedInterior( asQP( x ), u ); }

        virtual bool hasBoundary () const override { return impl().hasBoundary(); }
        virtual ValueType boundary ( const SurfaceCachingPointType &x, const ValueType &u ) const override { return impl().boundary( asQP( x ), u ); }
        virtual ValueType boundary ( const SurfaceElementPointType &x, const ValueType &u ) const override { return impl().boundary( asQP( x ), u ); }
        virtual Linearization< ValueType > linearizedBoundary ( const SurfaceCachingPointType &x, const ValueType &u ) const override { return impl().linearizedBoundary( asQP( x ), u ); }
        virtual Linearization< ValueType > linearizedBoundary ( const SurfaceElementPointType &x, const ValueType &u ) const override { return impl().linearizedBoundary( asQP( x ), u ); }

        virtual bool hasSkeleton () const override { return impl().hasSkeleton(); }
        virtual std::pair< ValueType, ValueType > skeleton ( const SurfaceCachingPointType &xIn, const ValueType &uIn, const SurfaceCachingPointType &xOut, const ValueType &uOut ) const override { return impl().skeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual std::pair< ValueType, ValueType > skeleton ( const SurfaceElementPointType &xIn, const ValueType &uIn, const SurfaceElementPointType &xOut, const ValueType &uOut ) const override { return impl().skeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual LinearizationPair< ValueType > linearizedSkeleton ( const SurfaceCachingPointType &xIn, const ValueType &uIn, const SurfaceCachingPointType &xOut, const ValueType &uOut ) const override { return impl().linearizedSkeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual LinearizationPair< ValueType > linearizedSkeleton ( const SurfaceElementPointType &xIn, const ValueType &uIn, const SurfaceElementPointType &xOut, const ValueType &uOut ) const override { return impl().linearizedSkeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }

      private:
        const auto &impl () const { return std::cref( impl_ ).get(); }
        auto &impl () { return std::ref( impl_ ).get(); }

        Impl impl_;
      };

    public:
      template< class Integrands, std::enable_if_t< IntegrandsTraits< std::decay_t< Integrands > >::isFull && !std::is_same< std::decay_t< Integrands >, This >::value, int > = 0 >
      explicit VirtualizedIntegrands ( Integrands integrands )
        : impl_( new Implementation< Integrands >( std::move( integrands ) ) )
      {}

      template< class Integrands, std::enable_if_t< !IntegrandsTraits< Integrands >::isFull, int > = 0 >
      explicit VirtualizedIntegrands ( Integrands integrands )
        : VirtualizedIntegrands( FullIntegrands< std::decay_t< Integrands > >( std::move( integrands ) ) )
      {}

      VirtualizedIntegrands ( const This &other ) : impl_( other ? other.impl_->clone() : nullptr ) {}
      VirtualizedIntegrands ( This && ) = default;

      VirtualizedIntegrands &operator= ( const This &other ) { impl_.reset( other ? other.impl_->clone() : nullptr ); }
      VirtualizedIntegrands &operator= ( This && ) = default;

      explicit operator bool () const { return static_cast< bool >( impl_ ); }

      bool init ( const EntityType &entity ) { return impl_->init( entity ); }
      bool init ( const IntersectionType &intersection ) { return impl_->init( intersection ); }

      bool hasInterior () const { return impl_->hasInterior(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      ValueType interior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->interior( InteriorCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      ValueType interior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->interior( InteriorElementPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedInterior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->linearizedInterior( InteriorCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedInterior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->linearizedInterior( InteriorElementPointType( x ), u );
      }

      bool hasBoundary () const { return impl_->hasBoundary(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      ValueType boundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->boundary( SurfaceCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      ValueType boundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->boundary( SurfaceElementPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedBoundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->linearizedBoundary( SurfaceCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedBoundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const ValueType &u ) const
      {
        return impl_->linearizedBoundary( SurfaceElementPointType( x ), u );
      }

      bool hasSkeleton () const { return impl_->hasSkeleton(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      std::pair< ValueType, ValueType > skeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const ValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const ValueType &uOut ) const
      {
        return impl_->skeleton( SurfaceCachingPointType( xIn ), uIn, SurfaceCachingPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      std::pair< ValueType, ValueType > skeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const ValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const ValueType &uOut ) const
      {
        return impl_->skeleton( SurfaceElementPointType( xIn ), uIn, SurfaceElementPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedSkeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const ValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const ValueType &uOut ) const
      {
        return impl_->linearizedSkeleton( SurfaceCachingPointType( xIn ), uIn, SurfaceCachingPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedSkeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const ValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const ValueType &uOut ) const
      {
        return impl_->linearizedSkeleton( SurfaceElementPointType( xIn ), uIn, SurfaceElementPointType( xOut ), uOut );
      }

    private:
      std::unique_ptr< Interface > impl_;
    };



    // DiffusionModelIntegrands
    // ------------------------

    template< class Model >
    struct DiffusionModelIntegrands
    {
      typedef Model ModelType;
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
      typedef Model ModelType;
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
          beta_ = penalty_ * intersection.geometry().volume() / inside.geometry().volume();
          return (model().hasNeumanBoundary() && model().init( inside ));
        }
        else if( intersection.neighbor() )
        {
          const auto volIn = intersection.inside().geometry().volume();
          const auto volOut = intersection.outside().geometry().volume();
          beta_ = penalty_ * intersection.geometry().volume() / std::min( volIn, volOut );
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
        return ValueType( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linDiffusiveFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return ValueType( source, dFlux );
          };
      }

      template< class Point >
      ValueType boundary ( const Point &x, const ValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return ValueType( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const ValueType &u ) const
      {
        return [ this, x, u ] ( const ValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return ValueType( alpha, 0 );
          };
      }

      template< class Point >
      std::pair< ValueType, ValueType > skeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const EntityType inside = intersection().inside();
        const EntityType outside = intersection().outside();

        const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        ValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        model().init( outside );
        JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uJump ), dFluxPrimeOut );
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), 0, dFluxOut );
        dFluxPrimeOut -= dFluxOut;
        dFluxOut = 0;
        model().diffusiveFlux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uOut ), dFluxOut );

        model().init( inside );
        JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uJump ), dFluxPrimeIn );
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), 0, dFluxIn );
        dFluxPrimeIn -= dFluxIn;
        dFluxIn = 0;
        model().diffusiveFlux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uIn ), dFluxIn );

        RangeType int0 = std::get< 0 >( uJump );
        int0 *= beta_;
        dFluxIn += dFluxOut;
        dFluxIn.usmv( -half, normal, int0 );

        dFluxPrimeIn *= -half;
        dFluxPrimeOut *= -half;

        return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const ValueType &uIn, const Point &xOut, const ValueType &uOut ) const
      {
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        ValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        auto intIn = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiIn ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          ValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) -= std::get< 0 >( phiIn );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, 0, std::get< 1 >( phiJump ), dFluxPrimeOut );

          model().init( inside );
          JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiJump ), dFluxPrimeIn );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), 0, xIn, std::get< 0 >( phiIn ), 0, dFluxIn );
          dFluxPrimeIn -= dFluxIn;
          dFluxIn = 0;
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uIn ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), dFluxIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxIn.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
        };

        auto intOut = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const ValueType &phiOut ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          ValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) = std::get< 0 >( phiOut );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiJump ), dFluxPrimeOut );
          model().linDiffusiveFlux( std::get< 0 >( uOut ), 0, xOut, std::get< 0 >( phiOut ), 0, dFluxOut );
          dFluxPrimeOut -= dFluxOut;
          dFluxOut = 0;
          model().linDiffusiveFlux( std::get< 0 >( uOut ), std::get< 1 >( uOut ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), dFluxOut );

          model().init( inside );
          JacobianRangeType dFluxPrimeIn( 0 );
          model().linDiffusiveFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, 0, std::get< 1 >( phiJump ), dFluxPrimeIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxOut.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( ValueType( -int0, dFluxPrimeIn ), ValueType( int0, dFluxPrimeOut ) );
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

#endif // #ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
