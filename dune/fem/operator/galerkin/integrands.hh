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
        static std::true_type interior ( const Integrands &, decltype( std::declval< const Integrands & >().interior( std::declval< const InteriorQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type interior ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasInterior() ), bool >::value, int > = 0 >
        static std::true_type hasInterior ( const Integrands & );

        static std::false_type hasInterior ( ... );


        template< class Integrands >
        static std::true_type boundary ( const Integrands &, decltype( std::declval< const Integrands & >().boundary( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type boundary ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasBoundary() ), bool >::value, int > = 0 >
        static std::true_type hasBoundary ( const Integrands & );

        static std::false_type hasBoundary ( ... );


        template< class Integrands >
        static std::true_type skeleton ( const Integrands &, decltype( std::declval< const Integrands & >().skeleton( std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >(), std::declval< const SurfaceQuadraturePointType< Integrands > & >(), std::declval< const DomainValueType< Integrands > & >() ) ) * = nullptr );

        static std::false_type skeleton ( ... );

        template< class Integrands, std::enable_if_t< std::is_same< decltype( std::declval< const Integrands & >().hasSkeleton() ), bool >::value, int > = 0 >
        static std::true_type hasSkeleton ( const Integrands & );

        static std::false_type hasSkeleton ( ... );

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
        : integrands_( std::forward< Args >( args )... )
      {}

      bool init ( const EntityType &entity ) { return integrands().init( entity ); }
      bool init ( const IntersectionType &intersection ) { return integrands().init( intersection ); }

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

      // Perhaps add boundary classification method
#if 0
      IdentificationClass classifyBoundary(Intersection)
      {
        return blah;
      }
#endif

    private:
      decltype( auto ) integrands () { return std::ref( integrands_ ).get(); }
      decltype( auto ) integrands () const { return std::ref( integrands_ ).get(); }

      Integrands integrands_;

    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_INTEGRANDS_HH
