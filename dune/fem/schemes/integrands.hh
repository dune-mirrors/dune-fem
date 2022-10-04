#ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
#define DUNE_FEM_SCHEMES_INTEGRANDS_HH

#include <cassert>

#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/fempy/quadrature/cachingpoint.hh>
#include <dune/fempy/quadrature/elementpoint.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/common/explicitfieldvector.hh>

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
        using DomainValueType = typename std::decay_t< decltype( std::ref( std::declval< const Integrands & >() ).get() ) >::DomainValueType;

        template< class Integrands >
        using RangeValueType = typename std::decay_t< decltype( std::ref( std::declval< const Integrands & >() ).get() ) >::RangeValueType;

        template< class Integrands >
        using GridPartType = typename std::decay_t< decltype( std::ref( std::declval< const Integrands & >() ).get() ) >::GridPartType;


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
        std::true_type interior ( const Integrands &, decltype( std::declval< const Integrands & >().interior( std::declval< const InteriorQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        std::false_type interior ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasInterior() ), bool >::value, int > = 0 >
        std::true_type hasInterior ( const Integrands & );

        std::false_type hasInterior ( ... );


        template< class Integrands >
        std::true_type boundary ( const Integrands &, decltype( std::declval< const Integrands & >().boundary( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        std::false_type boundary ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasBoundary() ), bool >::value, int > = 0 >
        std::true_type hasBoundary ( const Integrands & );

        std::false_type hasBoundary ( ... );


        template< class Integrands >
        std::true_type skeleton ( const Integrands &, decltype( std::declval< const Integrands & >().skeleton( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >(), std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        std::false_type skeleton ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasSkeleton() ), bool >::value, int > = 0 >
        std::true_type hasSkeleton ( const Integrands & );

        std::false_type hasSkeleton ( ... );

        template< class Integrands >
        using Get = std::decay_t< decltype( std::ref( std::declval< const Integrands & >() ).get() ) >;

        template< class Integrands >
        using RRangeType = typename Get<Integrands>::RRangeType;
        template< class Integrands >
        using DirichletComponentType = typename Get<Integrands>::DirichletComponentType;
      } // namespace IntegrandsTraits

    } // namespace Impl

    template< class Integrands >
    struct IntegrandsTraits
    {
      typedef Impl::IntegrandsTraits::DomainValueType< Integrands > DomainValueType;
      typedef Impl::IntegrandsTraits::RangeValueType< Integrands > RangeValueType;
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

      typedef Impl::IntegrandsTraits::RRangeType< Integrands > RRangeType;
      typedef Impl::IntegrandsTraits::DirichletComponentType< Integrands > DirichletComponentType;
    };



    // FullIntegrands
    // --------------

    template< class Integrands >
    struct FullIntegrands
    {
      typedef typename IntegrandsTraits< Integrands >::DomainValueType DomainValueType;
      typedef typename IntegrandsTraits< Integrands >::RangeValueType RangeValueType;
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
      static RangeValueType interior ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return integrands.interior( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::interior, int > = 0 >
      static RangeValueType interior ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return RangeValueType();
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::interior, int > = 0 >
      static auto linearizedInterior ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return integrands.linearizedInterior( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::interior, int > = 0 >
      static auto linearizedInterior ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return [] ( const DomainValueType & ) { return RangeValueType(); };
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
      static RangeValueType boundary ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return integrands.boundary( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::boundary, int > = 0 >
      static RangeValueType boundary ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return RangeValueType();
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::boundary, int > = 0 >
      static auto linearizedBoundary ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return integrands.linearizedBoundary( x, u );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::boundary, int > = 0 >
      static auto linearizedBoundary ( const T &integrands, const Point &x, const DomainValueType &u )
      {
        return [] ( const DomainValueType & ) { return RangeValueType(); };
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
      static std::pair< RangeValueType, RangeValueType > skeleton ( const T &integrands, const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut )
      {
        return integrands.skeleton( xIn, uIn, xOut, uOut );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::skeleton, int > = 0 >
      static std::pair< RangeValueType, RangeValueType > skeleton ( const T &integrands, const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut )
      {
        return std::make_pair( RangeValueType(), RangeValueType() );
      }

      template< class T, class Point, std::enable_if_t< IntegrandsTraits< T >::skeleton, int > = 0 >
      static auto linearizedSkeleton ( const T &integrands, const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut )
      {
        return integrands.linearizedSkeleton( xIn, uIn, xOut, uOut );
      }

      template< class T, class Point, std::enable_if_t< !IntegrandsTraits< T >::skeleton, int > = 0 >
      static auto linearizedSkeleton ( const T &integrands, const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut )
      {
        auto zero = [] ( const DomainValueType & ) { return std::make_pair( RangeValueType(), RangeValueType() ); };
        return std::make_pair( zero, zero );
      }

    public:
      template< class... Args >
      explicit FullIntegrands ( Args &&... args )
        : integrands_( std::forward< Args >( args )... ),
          rInt_( std::ref( integrands_ ).get() )
      {
      }

      bool init ( const EntityType &entity ) { return integrands().init( entity ); }
      bool init ( const IntersectionType &intersection ) { return integrands().init( intersection ); }
      void unbind ( ) { integrands().unbind( ); }

      bool hasInterior () const { return hasInterior( integrands() ); }

      template< class Point >
      RangeValueType interior ( const Point &x, const DomainValueType &u ) const
      {
        return interior( integrands(), x, u );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
      {
        return linearizedInterior( integrands(), x, u );
      }

      bool hasBoundary () const { return hasBoundary( integrands() ); }

      template< class Point >
      RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
      {
        return boundary( integrands(), x, u );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
      {
        return linearizedBoundary( integrands(), x, u );
      }

      bool hasSkeleton () const { return hasSkeleton( integrands() ); }

      template< class Point >
      std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
      {
        return skeleton( integrands(), xIn, uIn, xOut, uOut );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
      {
        return linearizedSkeleton( integrands(), xIn, uIn, xOut, uOut );
      }

    protected:
      typedef typename Integrands::type RealIntegrands;
      //decltype( auto ) integrands () { return std::ref( integrands_ ).get(); }
      //decltype( auto ) integrands () const { return std::ref( integrands_ ).get(); }
      RealIntegrands& integrands () { return rInt_; }
      const RealIntegrands& integrands () const { return rInt_; }



      Integrands integrands_;
      RealIntegrands rInt_;

    public:
      typedef typename IntegrandsTraits< Integrands >::RRangeType RRangeType;
      typedef typename IntegrandsTraits< Integrands >::DirichletComponentType DirichletComponentType;
      bool hasDirichletBoundary () const
      {
        return integrands().hasDirichletBoundary();
      }
      bool isDirichletIntersection( const IntersectionType& inter, DirichletComponentType &dirichletComponent ) const
      {
        return integrands().isDirichletIntersection(inter,dirichletComponent);
      }
      template <class Point>
      void dirichlet( int bndId, const Point &x, RRangeType &value) const
      {
        return integrands().dirichlet(bndId,x,value);
      }
    };



    // VirtualizedIntegrands
    // ---------------------

    namespace detail
    {
      template <class T>
      struct GetDimRange;
      template <class FT,int r>
      struct GetDimRange<Dune::FieldVector<FT,r>>
      {
        typedef Dune::FieldVector<FT,r> type;
        static const int value = r;
      };
      template <class FT,int r,int c>
      struct GetDimRange<Dune::FieldMatrix<FT,r,c>>
      {
        typedef Dune::FieldVector<FT,r> type;
        static const int value = r;
      };
      template <class FT,int r,int c>
      struct GetDimRange<Dune::Fem::ExplicitFieldVector<Dune::FieldMatrix<FT,c,c>,r>>
      {
        typedef Dune::FieldVector<FT,r> type;
        static const int value = r;
      };
    }

    template< class GridPart, class DomainValue, class RangeValue = DomainValue >
    class VirtualizedIntegrands
    {
      typedef VirtualizedIntegrands< GridPart, DomainValue, RangeValue > This;

    public:
      typedef GridPart GridPartType;
      typedef DomainValue DomainValueType;
      typedef RangeValue RangeValueType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      using RRangeType = typename detail::GetDimRange<std::tuple_element_t<0,RangeValueType>>::type;
      typedef std::array<int,RRangeType::dimension> DirichletComponentType;
      typedef typename EntityType::Geometry::LocalCoordinate DomainType;

    private:
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      typedef FemPy::CachingPoint< GridPart, LocalCoordinateType, 0 > InteriorCachingPointType;
      typedef FemPy::ElementPoint< GridPart, LocalCoordinateType, 0 > InteriorElementPointType;
      typedef FemPy::CachingPoint< GridPart, LocalCoordinateType, 1 > SurfaceCachingPointType;
      typedef FemPy::ElementPoint< GridPart, LocalCoordinateType, 1 > SurfaceElementPointType;

      template< class QP >
      static Fem::QuadraturePointWrapper< QP > asQP ( const QP &qp )
      {
        return static_cast< Fem::QuadraturePointWrapper< QP > >( qp );
      }

      template< class R >
      using Linearization = std::function< R( const DomainValueType & ) >;

      template< class R >
      using LinearizationPair = std::pair< Linearization< std::pair< R, R > >, Linearization< std::pair< R, R > > >;

      struct Interface
      {
        virtual ~Interface ()  {}
        virtual Interface *clone () const = 0;

        virtual bool init ( const EntityType &entity ) = 0;
        virtual bool init ( const IntersectionType &intersection ) = 0;
        virtual void unbind ( ) = 0;

        virtual bool hasInterior () const = 0;
        virtual RangeValueType interior ( const InteriorCachingPointType &x, const DomainValueType &u ) const = 0;
        virtual RangeValueType interior ( const InteriorElementPointType &x, const DomainValueType &u ) const = 0;
        virtual Linearization< RangeValueType > linearizedInterior ( const InteriorCachingPointType &x, const DomainValueType &u ) const = 0;
        virtual Linearization< RangeValueType > linearizedInterior ( const InteriorElementPointType &x, const DomainValueType &u ) const = 0;

        virtual bool hasBoundary () const = 0;
        virtual RangeValueType boundary ( const SurfaceCachingPointType &x, const DomainValueType &u ) const = 0;
        virtual RangeValueType boundary ( const SurfaceElementPointType &x, const DomainValueType &u ) const = 0;
        virtual Linearization< RangeValueType > linearizedBoundary ( const SurfaceCachingPointType &x, const DomainValueType &u ) const = 0;
        virtual Linearization< RangeValueType > linearizedBoundary ( const SurfaceElementPointType &x, const DomainValueType &u ) const = 0;

        virtual bool hasSkeleton () const = 0;
        virtual std::pair< RangeValueType, RangeValueType > skeleton ( const SurfaceCachingPointType &xIn, const DomainValueType &uIn, const SurfaceCachingPointType &xOut, const DomainValueType &uOut ) const = 0;
        virtual std::pair< RangeValueType, RangeValueType > skeleton ( const SurfaceElementPointType &xIn, const DomainValueType &uIn, const SurfaceElementPointType &xOut, const DomainValueType &uOut ) const = 0;
        virtual LinearizationPair< RangeValueType > linearizedSkeleton ( const SurfaceCachingPointType &xIn, const DomainValueType &uIn, const SurfaceCachingPointType &xOut, const DomainValueType &uOut ) const = 0;
        virtual LinearizationPair< RangeValueType > linearizedSkeleton ( const SurfaceElementPointType &xIn, const DomainValueType &uIn, const SurfaceElementPointType &xOut, const DomainValueType &uOut ) const = 0;

        virtual bool hasDirichletBoundary () const = 0;
        virtual bool isDirichletIntersection( const IntersectionType& inter, DirichletComponentType &dirichletComponent ) const = 0;
        virtual void dirichlet( int bndId, const DomainType &x,RRangeType &value) const = 0;
      };

      template< class Impl >
      struct DUNE_PRIVATE Implementation final
        : public Interface
      {
        Implementation ( Impl impl ) : impl_( std::move( impl ) )
        {
        }
        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual bool init ( const EntityType &entity ) override { return impl().init( entity ); }
        virtual bool init ( const IntersectionType &intersection ) override { return impl().init( intersection ); }
        virtual void unbind ( ) override { impl().unbind( ); }

        virtual bool hasInterior () const override { return impl().hasInterior(); }
        virtual RangeValueType interior ( const InteriorCachingPointType &x, const DomainValueType &u ) const override { return impl().interior( asQP( x ), u ); }
        virtual RangeValueType interior ( const InteriorElementPointType &x, const DomainValueType &u ) const override { return impl().interior( asQP( x ), u ); }
        virtual Linearization< RangeValueType > linearizedInterior ( const InteriorCachingPointType &x, const DomainValueType &u ) const override { return impl().linearizedInterior( asQP( x ), u ); }
        virtual Linearization< RangeValueType > linearizedInterior ( const InteriorElementPointType &x, const DomainValueType &u ) const override { return impl().linearizedInterior( asQP( x ), u ); }

        virtual bool hasBoundary () const override { return impl().hasBoundary(); }
        virtual RangeValueType boundary ( const SurfaceCachingPointType &x, const DomainValueType &u ) const override { return impl().boundary( asQP( x ), u ); }
        virtual RangeValueType boundary ( const SurfaceElementPointType &x, const DomainValueType &u ) const override { return impl().boundary( asQP( x ), u ); }
        virtual Linearization< RangeValueType > linearizedBoundary ( const SurfaceCachingPointType &x, const DomainValueType &u ) const override { return impl().linearizedBoundary( asQP( x ), u ); }
        virtual Linearization< RangeValueType > linearizedBoundary ( const SurfaceElementPointType &x, const DomainValueType &u ) const override { return impl().linearizedBoundary( asQP( x ), u ); }

        virtual bool hasSkeleton () const override { return impl().hasSkeleton(); }
        virtual std::pair< RangeValueType, RangeValueType > skeleton ( const SurfaceCachingPointType &xIn, const DomainValueType &uIn, const SurfaceCachingPointType &xOut, const DomainValueType &uOut ) const override { return impl().skeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual std::pair< RangeValueType, RangeValueType > skeleton ( const SurfaceElementPointType &xIn, const DomainValueType &uIn, const SurfaceElementPointType &xOut, const DomainValueType &uOut ) const override { return impl().skeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual LinearizationPair< RangeValueType > linearizedSkeleton ( const SurfaceCachingPointType &xIn, const DomainValueType &uIn, const SurfaceCachingPointType &xOut, const DomainValueType &uOut ) const override { return impl().linearizedSkeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }
        virtual LinearizationPair< RangeValueType > linearizedSkeleton ( const SurfaceElementPointType &xIn, const DomainValueType &uIn, const SurfaceElementPointType &xOut, const DomainValueType &uOut ) const override { return impl().linearizedSkeleton( asQP( xIn ), uIn, asQP( xOut ), uOut ); }

        virtual bool hasDirichletBoundary () const override { return impl().hasDirichletBoundary(); }
        virtual bool isDirichletIntersection( const IntersectionType& inter, DirichletComponentType &dirichletComponent ) const override { return impl().isDirichletIntersection(inter,dirichletComponent); }
        virtual void dirichlet( int bndId, const DomainType &x,RRangeType &value) const override { impl().dirichlet(bndId,x,value); }

      private:
        const auto &impl () const { return std::cref( impl_ ).get(); }
        auto &impl () { return std::ref( impl_ ).get(); }

        Impl impl_;
      };

      template< class Integrands >
      using isVirtualized = std::is_same< std::decay_t< decltype( std::ref( std::declval< Integrands & >() ).get() ) >, This >;

    public:
      template< class Integrands, std::enable_if_t< IntegrandsTraits< std::decay_t< Integrands > >::isFull && !isVirtualized< Integrands >::value, int > = 0 >
      explicit VirtualizedIntegrands ( Integrands integrands )
        : impl_( new Implementation< Integrands >( std::move( integrands ) ) )
      {}

      template< class Integrands, std::enable_if_t< !IntegrandsTraits< Integrands >::isFull && !isVirtualized< Integrands >::value, int > = 0 >
      explicit VirtualizedIntegrands ( Integrands integrands )
        : VirtualizedIntegrands( FullIntegrands< std::decay_t< Integrands > >( std::move( integrands ) ) )
      {}

      VirtualizedIntegrands ( const This &other ) : impl_( other ? other.impl().clone() : nullptr )
      {}

      VirtualizedIntegrands ( This && ) = default;

      VirtualizedIntegrands &operator= ( const This &other )
      {
        impl_.reset( other ? other.impl().clone() : nullptr );
      }

      VirtualizedIntegrands &operator= ( This && ) = default;

      explicit operator bool () const { return static_cast< bool >( impl_ ); }

      bool init ( const EntityType &entity ) { return impl().init( entity ); }
      bool init ( const IntersectionType &intersection ) { return impl().init( intersection ); }
      void unbind ( ) { impl().unbind( ); }

      bool hasInterior () const { return impl().hasInterior(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      RangeValueType interior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().interior( InteriorCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      RangeValueType interior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().interior( InteriorElementPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedInterior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().linearizedInterior( InteriorCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedInterior ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().linearizedInterior( InteriorElementPointType( x ), u );
      }

      bool hasBoundary () const { return impl().hasBoundary(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      RangeValueType boundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().boundary( SurfaceCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      RangeValueType boundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().boundary( SurfaceElementPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedBoundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().linearizedBoundary( SurfaceCachingPointType( x ), u );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedBoundary ( const Fem::QuadraturePointWrapper< Quadrature > &x, const DomainValueType &u ) const
      {
        return impl().linearizedBoundary( SurfaceElementPointType( x ), u );
      }

      bool hasSkeleton () const { return impl().hasSkeleton(); }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      std::pair< RangeValueType, RangeValueType > skeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const DomainValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const DomainValueType &uOut ) const
      {
        return impl().skeleton( SurfaceCachingPointType( xIn ), uIn, SurfaceCachingPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      std::pair< RangeValueType, RangeValueType > skeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const DomainValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const DomainValueType &uOut ) const
      {
        return impl().skeleton( SurfaceElementPointType( xIn ), uIn, SurfaceElementPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedSkeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const DomainValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const DomainValueType &uOut ) const
      {
        return impl().linearizedSkeleton( SurfaceCachingPointType( xIn ), uIn, SurfaceCachingPointType( xOut ), uOut );
      }

      template< class Quadrature, std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      auto linearizedSkeleton ( const Fem::QuadraturePointWrapper< Quadrature > &xIn, const DomainValueType &uIn, const Fem::QuadraturePointWrapper< Quadrature > &xOut, const DomainValueType &uOut ) const
      {
        return impl().linearizedSkeleton( SurfaceElementPointType( xIn ), uIn, SurfaceElementPointType( xOut ), uOut );
      }

      bool hasDirichletBoundary () const
      {
        return impl().hasDirichletBoundary();
      }
      bool isDirichletIntersection( const IntersectionType& inter, DirichletComponentType &dirichletComponent )  const
      {
        return impl().isDirichletIntersection(inter,dirichletComponent);
      }
      void dirichlet( int bndId, const DomainType &x,RRangeType &value) const
      {
        return impl().dirichlet(bndId,x,value);
      }

    private:
      const Interface &impl () const { assert( impl_ ); return *impl_; }
      Interface &impl () { assert( impl_ ); return *impl_; }

      std::unique_ptr< Interface > impl_;
    };



    // ConservationLawModelIntegrands
    // ------------------------------

    template< class Model >
    struct ConservationLawModelIntegrands
    {
      typedef Model ModelType;
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef std::tuple< RangeType, JacobianRangeType > DomainValueType;
      typedef std::tuple< RangeType, JacobianRangeType > RangeValueType;

      explicit ConservationLawModelIntegrands ( const Model &model ) : model_( &model ) {}

      bool init ( const EntityType &entity ) { return model().init( entity ); }

      bool init ( const IntersectionType &intersection )
      {
        return (intersection.boundary() && model().hasNeumanBoundary() && model().init( intersection.inside() ));
      }

      void unbind ( ) { model().unbind( ); }

      template< class Point >
      RangeValueType interior ( const Point &x, const DomainValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().flux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return std::make_tuple( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
      {
        return [ this, x, u ] ( const DomainValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return std::make_tuple( source, dFlux );
          };
      }

      template< class Point >
      RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return std::make_tuple( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
      {
        return [ this, x, u ] ( const DomainValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return std::make_tuple( alpha, 0 );
          };
      }

      const Model &model () const { return *model_; }

    private:
      const Model *model_ = nullptr;
    };



    // DGConservationLawModelIntegrands
    // --------------------------------

    template< class Model >
    struct DGConservationLawModelIntegrands
    {
      typedef Model ModelType;
      typedef typename Model::GridPartType GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename Model::RangeType RangeType;
      typedef typename Model::JacobianRangeType JacobianRangeType;

      typedef typename FieldTraits< RangeType >::field_type RangeFieldType;

      typedef std::tuple< RangeType, JacobianRangeType > DomainValueType;
      typedef std::tuple< RangeType, JacobianRangeType > RangeValueType;

      DGConservationLawModelIntegrands ( const Model &model, RangeFieldType penalty )
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

      void unbind ( )
      {
        model().unbind( );
      }

      template< class Point >
      RangeValueType interior ( const Point &x, const DomainValueType &u ) const
      {
        RangeType source( 0 );
        model().source( x, std::get< 0 >( u ), std::get< 1 >( u ), source );
        JacobianRangeType dFlux( 0 );
        model().flux( x, std::get< 0 >( u ), std::get< 1 >( u ), dFlux );
        return RangeValueType( source, dFlux );
      }

      template< class Point >
      auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
      {
        return [ this, x, u ] ( const DomainValueType &phi ) {
            RangeType source( 0 );
            model().linSource( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), source );
            JacobianRangeType dFlux( 0 );
            model().linFlux( std::get< 0 >( u ), std::get< 1 >( u ), x, std::get< 0 >( phi ), std::get< 1 >( phi ), dFlux );
            return RangeValueType( source, dFlux );
          };
      }

      template< class Point >
      RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
      {
        RangeType alpha( 0 );
        model().alpha( x, std::get< 0 >( u ), alpha );
        return RangeValueType( alpha, 0 );
      }

      template< class Point >
      auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
      {
        return [ this, x, u ] ( const DomainValueType &phi ) {
            RangeType alpha( 0 );
            model().linAlpha( std::get< 0 >( u ), x, std::get< 0 >( phi ), alpha );
            return RangeValueType( alpha, 0 );
          };
      }

      template< class Point >
      std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
      {
        const EntityType inside = intersection().inside();
        const EntityType outside = intersection().outside();

        const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        DomainValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        model().init( outside );
        JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
        model().flux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uJump ), dFluxPrimeOut );
        model().flux( xOut, std::get< 0 >( uOut ), 0, dFluxOut );
        dFluxPrimeOut -= dFluxOut;
        dFluxOut = 0;
        model().flux( xOut, std::get< 0 >( uOut ), std::get< 1 >( uOut ), dFluxOut );

        model().init( inside );
        JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
        model().flux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uJump ), dFluxPrimeIn );
        model().flux( xIn, std::get< 0 >( uIn ), 0, dFluxIn );
        dFluxPrimeIn -= dFluxIn;
        dFluxIn = 0;
        model().flux( xIn, std::get< 0 >( uIn ), std::get< 1 >( uIn ), dFluxIn );

        RangeType int0 = std::get< 0 >( uJump );
        int0 *= beta_;
        dFluxIn += dFluxOut;
        dFluxIn.usmv( -half, normal, int0 );

        dFluxPrimeIn *= -half;
        dFluxPrimeOut *= -half;

        return std::make_pair( RangeValueType( -int0, dFluxPrimeIn ), RangeValueType( int0, dFluxPrimeOut ) );
      }

      template< class Point >
      auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
      {
        const auto normal = intersection().unitOuterNormal( xIn.localPosition() );

        DomainValueType uJump( 0, 0 );
        std::get< 0 >( uJump ) = std::get< 0 >( uOut ) - std::get< 0 >( uIn );
        for( int i = 0; i < RangeType::dimension; ++i )
          std::get< 1 >( uJump )[ i ].axpy( std::get< 0 >( uJump )[ i ], normal );

        auto intIn = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const DomainValueType &phiIn ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          DomainValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) -= std::get< 0 >( phiIn );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxPrimeOut( 0 );
          model().linFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, 0, std::get< 1 >( phiJump ), dFluxPrimeOut );

          model().init( inside );
          JacobianRangeType dFluxIn( 0 ), dFluxPrimeIn( 0 );
          model().linFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiJump ), dFluxPrimeIn );
          model().linFlux( std::get< 0 >( uIn ), 0, xIn, std::get< 0 >( phiIn ), 0, dFluxIn );
          dFluxPrimeIn -= dFluxIn;
          dFluxIn = 0;
          model().linFlux( std::get< 0 >( uIn ), std::get< 1 >( uIn ), xIn, std::get< 0 >( phiIn ), std::get< 1 >( phiIn ), dFluxIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxIn.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( RangeValueType( -int0, dFluxPrimeIn ), RangeValueType( int0, dFluxPrimeOut ) );
        };

        auto intOut = [ this, xIn, uIn, xOut, uOut, normal, uJump ] ( const DomainValueType &phiOut ) {
          const EntityType inside = intersection().inside();
          const EntityType outside = intersection().outside();

          const RangeFieldType half = RangeFieldType( 1 ) / RangeFieldType( 2 );

          DomainValueType phiJump( 0, 0 );
          std::get< 0 >( phiJump ) = std::get< 0 >( phiOut );
          for( int i = 0; i < RangeType::dimension; ++i )
            std::get< 1 >( phiJump )[ i ].axpy( std::get< 0 >( phiJump )[ i ], normal );

          model().init( outside );
          JacobianRangeType dFluxOut( 0 ), dFluxPrimeOut( 0 );
          model().linFlux( std::get< 0 >( uOut ), std::get< 1 >( uJump ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiJump ), dFluxPrimeOut );
          model().linFlux( std::get< 0 >( uOut ), 0, xOut, std::get< 0 >( phiOut ), 0, dFluxOut );
          dFluxPrimeOut -= dFluxOut;
          dFluxOut = 0;
          model().linFlux( std::get< 0 >( uOut ), std::get< 1 >( uOut ), xOut, std::get< 0 >( phiOut ), std::get< 1 >( phiOut ), dFluxOut );

          model().init( inside );
          JacobianRangeType dFluxPrimeIn( 0 );
          model().linFlux( std::get< 0 >( uIn ), std::get< 1 >( uJump ), xIn, 0, std::get< 1 >( phiJump ), dFluxPrimeIn );

          RangeType int0 = std::get< 0 >( phiJump );
          int0 *= beta_;
          dFluxOut.usmv( -half, normal, int0 );

          dFluxPrimeIn *= -half;
          dFluxPrimeOut *= -half;

          return std::make_pair( RangeValueType( -int0, dFluxPrimeIn ), RangeValueType( int0, dFluxPrimeOut ) );
        };

        return std::make_pair( intIn, intOut );
      }

      const Model &model () const { return *model_; }

    private:
      const IntersectionType &intersection () const { assert( intersection_ ); return *intersection_; }

      const Model *model_ = nullptr;
      RangeFieldType penalty_;
      const IntersectionType *intersection_ = nullptr;
      RangeFieldType beta_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
