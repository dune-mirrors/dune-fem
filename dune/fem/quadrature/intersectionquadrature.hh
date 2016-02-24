#ifndef DUNE_FEM_INTERSECTIONQUADRATURE_HH
#define DUNE_FEM_INTERSECTIONQUADRATURE_HH

//- Dune includes
#include <dune/common/math.hh>

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

#include "elementquadrature.hh"
#include "cachingquadrature.hh"

namespace Dune
{

  namespace Fem
  {

    /** \brief IntersectionQuadrature is a helper class for creating the appropriate face quadratures
               for integrating over intersections. */
    template< typename FaceQuadrature, bool conforming  >
    class IntersectionQuadrature
    {
      template < typename FaceQuadratureImp, bool isConforming >
      struct QuadSelector
      {
        // use given quadrature
        typedef FaceQuadratureImp  FaceQuadratureType;
      };

      template < typename FaceQuadratureImp >
      struct QuadSelector<FaceQuadratureImp, false>
      {
        // in this case non conforming type is used
        typedef typename  FaceQuadratureImp ::
          NonConformingQuadratureType  FaceQuadratureType;
      };

    public:
      //! type of grid partition
      typedef typename FaceQuadrature :: GridPartType GridPartType;

      //! Type of the intersection iterator
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      //! codimension of the element quadrature
      enum { codimension = FaceQuadrature :: codimension };

      //! type of intersection quadrature implementation
      typedef typename QuadSelector<FaceQuadrature, conforming> :: FaceQuadratureType FaceQuadratureType;

      //! Dimension of the world.
      enum { dimension = FaceQuadratureType ::dimension };

      //! Just another name for double...
      typedef typename FaceQuadratureType :: RealType RealType;
      //! The type of the coordinates in the codim-0 reference element.
      typedef typename FaceQuadratureType :: CoordinateType CoordinateType;

      typedef typename FaceQuadratureType::LocalCoordinateType LocalCoordinateType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      /** \brief Constructor creating an inside and an outside face quadrature for
                 integrating over an intersection.

          \param[in]  gridPart      grid partition
          \param[in]  intersection  intersection
          \param[in]  order         desired order of the quadrature

       */
      IntersectionQuadrature( const GridPartType &gridPart,
                              const IntersectionType &intersection,
                              const int order)
      : inside_ ( gridPart, intersection, order, FaceQuadratureType::INSIDE ),
        outside_( gridPart, intersection, order, intersection.neighbor() ? FaceQuadratureType::OUTSIDE : FaceQuadratureType::INSIDE )
      {}

      /** \brief Constructor creating an inside and an outside face quadrature for
                 integrating over an intersection.

          \param[in]  gridPart        grid partition
          \param[in]  intersection    intersection
          \param[in]  order           desired order of the quadrature
          \param[in]  noNeigborCheck  flag that indicates that the neighbor check is not necessary (independent of the value of noNeigborCheck)

          \note For this constructor intersection.neighbor() must return true.
       */
      IntersectionQuadrature( const GridPartType &gridPart,
                              const IntersectionType &intersection,
                              const int order,
                              const bool noNeighborCheck )
      : inside_ ( gridPart, intersection, order, FaceQuadratureType::INSIDE  ),
        outside_( gridPart, intersection, order, FaceQuadratureType::OUTSIDE )
      {
        // make sure neighbor is true
        assert( intersection.neighbor() );
      }

      //! \brief return reference to inside face quadrature
      const FaceQuadratureType& inside()  const { return inside_;  }

      //! \brief return reference to outside face quadrature
      const FaceQuadratureType& outside() const { return outside_; }

      size_t nop () const
      {
        assert( inside().nop() == outside().nop() );
        return inside().nop();
      }

      const LocalCoordinateType &localPoint ( const int qp ) const
      {
        assert( inside().localPoint( qp ) == outside().localPoint( qp ) );
        return inside().localPoint( qp );
      }

      const RealType &weight ( const int qp ) const
      {
        assert( inside().weight( qp ) == outside().weight( qp ) );
        return inside().weight( qp );
      }

      IntersectionQuadrature( const IntersectionQuadrature& ) = delete;

    protected:
      const FaceQuadratureType inside_;
      const FaceQuadratureType outside_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_INTERSECTIONQUADRATURE_HH
