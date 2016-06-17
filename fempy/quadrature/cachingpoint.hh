#ifndef DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH
#define DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH

#include <type_traits>

#include <dune/common/ftraits.hh>

#include <dune/fem/quadrature/cachingpointlist.hh>

namespace Dune
{

  namespace FemPy
  {

    // CachingPoint
    // ------------

    template< class Coordinate >
    struct CachingPoint
      : public Fem::CachingInterface
    {
      typedef CachingPoint< Coordinate > This;

    public:
      typedef Coordinate CoordinateType;
      typedef typename FieldTraits< Coordinate >::real_type RealType;
      typedef Coordinate LocalCoordinateType;

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      CachingPoint ( const Quadrature &quadrature, std::size_t idx )
        : id_( quadrature.id() ),
          idx_( quadrature.cachingPoint( idx ) ),
          position_( quadrature.point( idx ) )
      {}

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      explicit CachingPoint ( const Fem::QuadraturePointWrapper< Quadrature > &x )
#if DUNE_VERSION_NEWER_REV(DUNE_FEM,2,5,0)
        : CachingPoint( x.quadrature(), x.index() )
#else
        : CachingPoint( x.quadrature(), x.point() )
#endif
      {}

      explicit operator Fem::QuadraturePointWrapper< This > () const noexcept { return Fem::QuadraturePointWrapper< This >( *this, 0u ); }

      std::size_t id () const { return id_; }

      const CoordinateType &point ( std::size_t qp ) const { return position_; }

      std::size_t cachingPoint ( std::size_t qp ) const { assert( qp == 0u ); return idx_; }

    private:
      std::size_t id_, idx_;
      const Coordinate &position_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH
