#ifndef DUNE_FEM_QUADRATURE_LIFTED_HH
#define DUNE_FEM_QUADRATURE_LIFTED_HH

#include <cassert>

#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/quadrature.hh>

namespace Dune
{

  namespace Fem
  {

    // IntersectionSide
    // ----------------

    enum class IntersectionSide { inside, outside };



    // FakeGridPart
    // ------------

    template< int dim >
    struct FakeGridPart
    {
      typedef double ctype;
      static const int dimension = dim;

      template< int codim >
      struct Codim
      {
        typedef int EntityType;
      };
    };



    // LiftedQuadrature
    // ----------------

    template< class Quadrature >
    class LiftedQuadrature
    {
      typedef LiftedQuadrature< Quadrature > ThisType;

    public:
      typedef typename Quadrature::CoordinateType CoordinateType;
      typedef typename Quadrature::LocalCoordinateType LocalCoordinateType;

      typedef QuadratureIterator< ThisType > IteratorType;
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;

      template< class Lifting >
      LiftedQuadrature ( Quadrature quadrature, const Lifting &lifting )
        : quadrature_( std::move( quadrature ) ), weights_( quadrature.nop() )
      {
        for( std::size_t i = 0; i < weights_.size(); ++i )
          weights_[ i ] = lifting( quadrature[ i ] );
      }

      QuadraturePointWrapperType operator[] ( std::size_t i ) const noexcept { return QuadraturePointWrapperType( *this, i ); }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      std::size_t cachingPoint ( std::size_t i ) const { return base().cachingPoint( i ); }
      const LocalCoordinateType &localPoint ( std::size_t i ) const { return base().localPoint( i ); }
      const CoordinateType &point ( std::size_t i ) const { return base().point( i ); }

      std::size_t id () const { return base().id(); }
      std::size_t nop () const { return base().nop(); }
      int order () const { return base().order(); }
      GeometryType type () const { return base().type(); }

      const RealType &weight ( std::size_t i ) const { return weights_[ i ]; }

      const Quadrature &base () const noexcept { return quadrature_; }

    private:
      Quadrature quadrature_;
      std::vector< typename Quadrature::RealType > weights_;
    };



    // liftedQuadrature
    // ----------------

    template< class Quadrature, class Geometry >
    inline static auto liftedQuadrature ( Quadrature quadrature, const Geometry &geometry )
    {
      if( geometry.affine() )
      {
        const auto integrationElement = geometry.integrationElement( CoordinateType( 0 ) );
        auto lifting = [ integrationElement ] ( const QuadraturePointWrapper< Quadrature > &x ) {
            return x.weight() * integrationElement;
          };
        return LiftingQuadrature< Quadrature >( std::move( quadrature ), lifting );
      }
      else
      {
        auto lifting = [ &geometry ] ( const QuadraturePointWrapper< Quadrature > &x ) {
            return x.weight() * geometry.integrationElement( x.position() );
          };
        return LiftingQuadrature< Quadrature >( std::move( quadrature ), lifting );
      }
    }



    // liftedCachingQuadrature
    // -----------------------

    template< int dim, int codim, class Grid, class Impl, class Lifting >
    inline static auto liftedCachingQuadrature ( const Entity< dim, codim, Grid, Impl > &entity, int order, const Lifting &lifting )
    {
      typedef CachingQuadrature< FakeGridPart< dim >, codim > Quadrature;
      return LiftingQuadrature< Quadrature >( Quadrature( entity, order ), lifting );
    }

    template< int dim, int codim, class Grid, class Impl >
    inline static auto liftedCachingQuadrature ( const Entity< dim, codim, Grid, Impl > &entity, int order )
    {
      typedef CachingQuadrature< FakeGridPart< dim >, codim > Quadrature;
      return liftedQuadrature( Quadrature( entity, order ), entity.geometry() );
    }

    template< class GridPart, class Lifting >
    inline static auto liftedCachingQuadrature ( const GridPart &gridPart, const typename GridPart::IntersectionType &intersection, int order, IntersectionSide side, const Lifting &lifting )
    {
      typedef CachingQuadrature< GridPart, 1 > Quadrature;
      switch( side )
      {
      case IntersectionSide::inside:
        return LiftedQuadrature< Quadrature >( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), lifting );
      case IntersectionSide::outside:
        return LiftedQuadrature< Quadrature >( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), lifting );
    }

    template< class GridPart, class Lifting >
    inline static auto liftedCachingQuadrature ( const GridPart &gridPart, const typename GridPart::IntersectionType &intersection, int order, IntersectionSide side )
    {
      typedef CachingQuadrature< GridPart, 1 > Quadrature;
      switch( side )
      {
      case IntersectionSide::inside:
        return liftedQuadrature( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), intersection.geometry() );
      case IntersectionSide::outside:
        return liftedQuadrature( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), intersection.geometry() );
    }



    // liftedElementQuadrature
    // -----------------------

    template< int dim, int codim, class Grid, class Impl, class Lifting >
    inline static auto liftedElementQuadrature ( const Entity< dim, codim, Grid, Impl > &entity, int order, const Lifting &lifting )
    {
      typedef ElementQuadrature< FakeGridPart< dim >, codim > Quadrature;
      return LiftingQuadrature< Quadrature >( Quadrature( entity, order ), lifting );
    }

    template< int dim, int codim, class Grid, class Impl >
    inline static auto liftedElementQuadrature ( const Entity< dim, codim, Grid, Impl > &entity, int order )
    {
      typedef ElementQuadrature< FakeGridPart< dim >, codim > Quadrature;
      return liftedQuadrature( Quadrature( entity, order ), entity.geometry() );
    }

    template< class GridPart, class Lifting >
    inline static auto liftedElementQuadrature ( const GridPart &gridPart, const typename GridPart::IntersectionType &intersection, int order, IntersectionSide side, const Lifting &lifting )
    {
      typedef ElementQuadrature< GridPart, 1 > Quadrature;
      switch( side )
      {
      case IntersectionSide::inside:
        return LiftedQuadrature< Quadrature >( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), lifting );
      case IntersectionSide::outside:
        return LiftedQuadrature< Quadrature >( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), lifting );
    }

    template< class GridPart, class Lifting >
    inline static auto liftedElementQuadrature ( const GridPart &gridPart, const typename GridPart::IntersectionType &intersection, int order, IntersectionSide side )
    {
      typedef ElementQuadrature< GridPart, 1 > Quadrature;
      switch( side )
      {
      case IntersectionSide::inside:
        return liftedQuadrature( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), intersection.geometry() );
      case IntersectionSide::outside:
        return liftedQuadrature( Quadrature( gridPart, intersection, order, Quadrature::OUTSIDE ), intersection.geometry() );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_QUADRATURE_LIFTED_HH
