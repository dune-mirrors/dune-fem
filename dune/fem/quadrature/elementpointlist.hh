#ifndef DUNE_FEM_ELEMENTPOINTLIST_HH
#define DUNE_FEM_ELEMENTPOINTLIST_HH

#include <dune/fem/quadrature/elementpointlistbase.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
{

  namespace Fem
  {

    /*! \class ElementIntegrationPointList
     *  \ingroup Quadrature
     *  \brief integration point list on the codim-0 reference element
     *
     *  DUNE quadratures are defined per geometry type, using local coordinates
     *  for the quadrature points. To evaluate a base function in some quadrature
     *  point, the quadrature must return points within the codim-0 reference
     *  element.
     *
     *  Now, assume you want to integrate over the face of a tetrahedron. This
     *  means you need a quadrature for a triangle, but the quadrature points
     *  should be specified with respect to the tetrahedron, since we want to
     *  evaluate our function in these points. This is where the ElementQuadrature
     *  comes into play.
     *
     *  The ElementIntegrationPointList takes a subentity and transforms the
     *  integration point list corresponding to the geometry to the codim-0
     *  reference element.
     *
     *  To achieve this goal, an ElementIntegrationPointList depends stronger on
     *  the context in which it is used. For example, for each face within a
     *  tetrahedron (though they are all the same) we need a different
     *  ElementIntegrationPointList, since the coordinates of the quadrature points
     *  with respect to the codim-0 entity differ for each face.
     *
     *  \note Actually, codim-1 element integration point lists depend on the
     *        intersection.
     *
     *  \note This integration point list does not support caching of base functions
     *        in integration points (see also CachingPointList).
     *
     *  For the actual implementations see
     *  - ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>
     *  - ElementIntegrationPointList<GridPartImp,1,IntegrationTraits>
     */
    template< class GridPartImp, int codim, class IntegrationTraits >
    class ElementIntegrationPointList;



    /** \copydoc ElementIntegrationPointList */
    template< class GridPartImp, class IntegrationTraits >
    class ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits >
    : public ElementPointListBase< GridPartImp, 0, IntegrationTraits >
    {
      typedef ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits > This;
      typedef ElementPointListBase< GridPartImp, 0, IntegrationTraits >  Base;

    public:
      //! type for coordinates in the codim-0 reference element
      typedef typename Base::CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename Base::QuadratureKeyType  QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< This > IteratorType;

      typedef typename Base :: IntegrationPointListType  IntegrationPointListType;

    public:
      using Base::localPoint;
      using Base::nop;

      /** \brief constructor
       *
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       */
      ElementIntegrationPointList( const GeometryType &geometry, const QuadratureKeyType& quadKey )
      : Base( geometry, quadKey )
      {}

      ElementIntegrationPointList( const IntegrationPointListType& ipList )
      : Base( ipList )
      {}

      const QuadraturePointWrapperType operator[] ( const size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /** \copydoc Dune::Fem::IntegrationPointList::point */
      const CoordinateType &point ( const size_t i ) const
      {
        return localPoint( i );
      }
    };



    /** \copydoc ElementIntegrationPointList */
    template< class GridPartImp, class IntegrationTraits >
    class ElementIntegrationPointList< GridPartImp, 1, IntegrationTraits >
    : public ElementPointListBase< GridPartImp, 1, IntegrationTraits >
    {
      typedef ElementIntegrationPointList< GridPartImp, 1, IntegrationTraits > This;
      typedef ElementPointListBase< GridPartImp, 1, IntegrationTraits > Base;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      static const int dimension = Base::dimension;

      //! Type of coordinates in codim-0 reference element
      typedef typename Base::CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename Base::QuadratureKeyType  QuadratureKeyType;

      //! Type of the intersection iterator
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< This > IteratorType;

      //! type quadrature for use on non-conforming intersections
      typedef This NonConformingQuadratureType;


      // for compatibility
      typedef typename GridPartType::TwistUtilityType  TwistUtilityType;
      typedef IntersectionIteratorType IntersectionIterator;


      using Base::localPoint;
      using Base::elementGeometry;
      using Base::nop;

      /** \brief constructor
       *
       *  \param[in]  gridPart      grid partition (a dummy here)
       *  \param[in]  intersection  intersection
       *  \param[in]  quadKey       quadrature key, i.e. desired order of the quadrature
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for
       *                            which the ElementQuadrature shall be created
       *
       *  \note This code assumes that the codim-0 entity is either a simplex or
       *        a cube (otherwise elementGeometry() returns a wrong geometry).
       */
      ElementIntegrationPointList ( const GridPartType &gridPart,
                                    const IntersectionType &intersection,
                                    const QuadratureKeyType& quadKey,
                                    const typename Base :: Side side )
      : Base( getPointList( intersection, quadKey, side ) ),
        referenceGeometry_( side == Base::INSIDE ?  intersection.geometryInInside() : intersection.geometryInOutside()),
        intersection_(intersection)
      {}

      const QuadraturePointWrapperType operator[] ( size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /** \copydoc Dune::Fem::IntegrationPointList::point
       */
      const CoordinateType &point ( size_t i ) const
      {
        dummy_ = referenceGeometry_.global( localPoint( i ) );
        return dummy_;
      }

      using Base::localFaceIndex;

      const IntersectionType &intersection() const
      {
        return intersection_;
      }
    protected:
      Base getPointList ( const IntersectionType &intersection, const int order,
                          const typename Base :: Side side )
      {
        switch( side )
        {
          case Base :: INSIDE:
            return Base( TwistUtilityType::elementGeometry( intersection, true ),
                         intersection.type(), intersection.indexInInside(), order );

          case Base ::OUTSIDE:
            return Base( TwistUtilityType::elementGeometry( intersection, false ),
                         intersection.type(), intersection.indexInOutside(), order );

          default:
            DUNE_THROW( InvalidStateException, "ElementIntegrationPointList: side must either be INSIDE or OUTSIDE." );
        }
      }

    private:
      typedef typename IntersectionIteratorType::Intersection::LocalGeometry ReferenceGeometry;

      ReferenceGeometry referenceGeometry_;
      const IntersectionType &intersection_;
      mutable CoordinateType dummy_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTPOINTLIST_HH
