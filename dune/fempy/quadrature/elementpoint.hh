#ifndef DUNE_FEMPY_QUADRATURE_ELEMENTPOINT_HH
#define DUNE_FEMPY_QUADRATURE_ELEMENTPOINT_HH

#include <type_traits>

#include <dune/common/ftraits.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune
{

  namespace FemPy
  {

    // ElementPoint
    // ------------

    template< class GridPart, class Coordinate, int codim >
    struct ElementPoint;

    template< class GridPart, class Coordinate >
    struct ElementPoint< GridPart, Coordinate, 0 >
    {
      typedef ElementPoint< GridPart, Coordinate, 0 > This;

    public:
      static const int codimension = 0 ;
      typedef Coordinate CoordinateType;
      typedef typename FieldTraits< Coordinate >::real_type RealType;
      typedef Coordinate LocalCoordinateType;

      template< class Quadrature >
      ElementPoint ( const Quadrature &quadrature, std::size_t idx )
        : position_( quadrature.point( idx ) )
      {}

      template< class Quadrature >
      explicit ElementPoint ( const Fem::QuadraturePointWrapper< Quadrature > &x )
        : ElementPoint( x.quadrature(), x.index() )
      {}

      explicit operator Fem::QuadraturePointWrapper< This > () const noexcept { return Fem::QuadraturePointWrapper< This >( *this, 0u ); }

      const CoordinateType &point ( std::size_t qp ) const { return position_; }
      const LocalCoordinateType &localPoint ( std::size_t qp ) const { return position_; }

    private:
      const CoordinateType &position_;
    };

    template< class GridPart, class Coordinate >
    struct ElementPoint< GridPart, Coordinate, 1 >
    {
      typedef ElementPoint< GridPart, Coordinate, 1 > This;

    public:
      static const int codimension = 1 ;
      typedef Coordinate CoordinateType;
      typedef typename FieldTraits< Coordinate >::real_type RealType;
      typedef FieldVector< typename FieldTraits< Coordinate >::field_type, Coordinate::dimension-1 > LocalCoordinateType;
      typedef typename GridPart::IntersectionType IntersectionType;

      template< class Quadrature >
      ElementPoint ( const Quadrature &quadrature, std::size_t idx )
        : position_( quadrature.point( idx ) ),
          localPosition_( quadrature.localPoint( idx ) ),
          intersection_( quadrature.intersection() )
      {}

      template< class Quadrature >
      explicit ElementPoint ( const Fem::QuadraturePointWrapper< Quadrature > &x )
        : ElementPoint( x.quadrature(), x.index() )
      {}

      explicit operator Fem::QuadraturePointWrapper< This > () const noexcept { return Fem::QuadraturePointWrapper< This >( *this, 0u ); }

      const CoordinateType &point ( std::size_t qp ) const { return position_; }
      const LocalCoordinateType &localPoint ( std::size_t qp ) const { return localPosition_; }

      const IntersectionType &intersection() const
      {
        return intersection_;
      }

    private:
      CoordinateType position_;
      const LocalCoordinateType &localPosition_;
      const IntersectionType &intersection_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_QUADRATURE_ELEMENTPOINT_HH
