#ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
#define DUNE_FEM_ELEMENTQUADRATURE_HH

#include <dune/fem/quadrature/elementpointlistbase.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
{

  namespace Fem
  {

    template< typename GridPartImp, class IntegrationPointList >
    class Agglomeration;

    template< class GridPartImp, int codim, template< class, int > class QuadratureTraits >
    struct ElementQuadratureTraits
    {
      // type of single coordinate
      typedef typename GridPartImp :: ctype ctype;

      // dimension of quadrature
      enum { dimension = GridPartImp ::  dimension };

      // codimension of quadrature
      enum { codimension = codim };

      // type of used integration point list
      typedef Quadrature< ctype, dimension-codim, QuadratureTraits > IntegrationPointListType;

      // type of local coordinate (with respect to the codim-0 entity)
      typedef typename Quadrature< ctype, dimension, QuadratureTraits > :: CoordinateType
        CoordinateType;
    };

    /*! \class ElementQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature on the codim-0 reference element
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
     *  The ElementQuadrature takes a subentity and transforms the quadrature
     *  corresponding to the geometry to the codim-0 reference element.
     *
     *  To achieve this goal, an element quadrature depends stronger on the
     *  context in which it is used. For example, for each face within a
     *  tetrahedron (though they are all the same) we need a different
     *  ElementQuadrature, since the coordinates of the quadrature points
     *  with respect to the codim-0 entity differ for each face.
     *
     *  \note Actually, codim-1 element quadratures depend on the intersection.
     *
     *  \note QuadratureTraits specifies the quadrature points used to build the
     *        quadrature.
     *
     *  \note This quadrature does not support caching of base functions in
     *        quadrature points (see also CachingQuadrature).
     *
     *  For the actual implementations, see
     *  - ElementQuadratureImpl<GridPartImp,0>
     *  - ElementQuadratureImpl<GridPartImp,1>
     */
    template< class GridPartImp, int codim, class IntegrationTraits, bool isQuadrature >
    class ElementQuadratureImpl;


    /** \copydoc ElementQuadratureImpl */
    template< class GridPartImp, class IntegrationTraits, bool isQuadrature >
    class ElementQuadratureImpl< GridPartImp, 0, IntegrationTraits, isQuadrature >
    : public ElementPointListBase< GridPartImp, 0, IntegrationTraits >
    {
      typedef ElementQuadratureImpl< GridPartImp, 0, IntegrationTraits, isQuadrature > This;
      typedef ElementPointListBase< GridPartImp, 0, IntegrationTraits >  Base;

    public:
      //! type of grid part
      typedef typename Base :: GridPartType  GridPartType;

      //! type for coordinates in the codim-0 reference element
      typedef typename Base::CoordinateType  CoordinateType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< This > IteratorType;

      // for compatibility
      typedef typename Base :: EntityType     EntityType;

      typedef typename Base :: IntegrationPointListType  IntegrationPointListType;

    protected:
      template <class QuadratureKeyType>
      IntegrationPointListType createQuadrature(  const EntityType &entity, const QuadratureKeyType& quadKey, const bool checkGeomType )
      {
        const GeometryType geomType = entity.type();
        if( checkGeomType && ! geomType.isNone() )
        {
          // return default element quadratures for given geometry type
          return IntegrationPointListType( geomType, quadKey );
        }
        else // compute weights and points based on sub-triangulation
        {
          typedef Agglomeration< GridPartType, IntegrationPointListType > AgglomerationType;
          return AgglomerationType::computeQuadrature( entity, quadKey );
        }
      }

      using Base::quadImp;

    public:
      using Base::localPoint;
      using Base::nop;

      /** \brief constructor
       *
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       *  \param[in]  checkGeomType if true geometry type is checked for isNone.
       */
      template <class QuadratureKeyType>
      ElementQuadratureImpl( const EntityType &entity, const QuadratureKeyType& quadKey, const bool checkGeomType = isQuadrature )
      : Base( createQuadrature( entity, quadKey, checkGeomType ) )
      {}

      /** \brief constructor
       *
       *  \param[in]  geometry  geometry type, the quadrature lives on
       *  \param[in]  order     desired minimal order of the quadrature
       */
      template <class QuadratureKeyType>
      ElementQuadratureImpl( const GeometryType &geometry, const QuadratureKeyType& quadKey )
      : Base( geometry, quadKey )
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

      /** \copydoc Dune::Fem::IntegrationPointList::weight */
      auto weight( size_t i ) const
      {
        return quadImp().weight( i );
      }
    };



    /** \copydoc ElementQuadratureImpl  */
    template< class GridPartImp, class IntegrationTraits, bool isQuadrature >
    class ElementQuadratureImpl< GridPartImp, 1, IntegrationTraits, isQuadrature >
    : public ElementPointListBase< GridPartImp, 1, IntegrationTraits >
    {
      typedef ElementQuadratureImpl< GridPartImp, 1, IntegrationTraits, isQuadrature > This;
      typedef ElementPointListBase< GridPartImp, 1, IntegrationTraits > Base;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      static const int dimension = Base::dimension;

      //! Type of coordinates in codim-0 reference element
      typedef typename Base::CoordinateType CoordinateType;

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
      template <class QuadratureKeyType>
      ElementQuadratureImpl( const GridPartType &gridPart,
                                       const IntersectionType &intersection,
                                       const QuadratureKeyType& quadKey,
                                       const typename Base :: Side side )
      : Base( getPointList( intersection, quadKey, side ) ),
        side_(side),
        intersection_(intersection),
        referenceGeometry_( side == Base::INSIDE ?  intersection.geometryInInside() : intersection.geometryInOutside())
      {}

      const QuadraturePointWrapperType operator[] ( size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }
      typename Base :: Side side() const { return side_; }
      bool isInside() const { return side_ == Base::INSIDE; }

      /** \copydoc Dune::Fem::IntegrationPointList::point */
      const CoordinateType &point ( size_t i ) const
      {
        dummy_ = referenceGeometry_.global( localPoint( i ) );
        return dummy_;
      }

      /** \copydoc Dune::Fem::IntegrationPointList::weight */
      auto weight( size_t i ) const
      {
        return quadImp().weight( i );
      }

      using Base::localFaceIndex;

      const IntersectionType &intersection() const
      {
        return intersection_;
      }

    protected:
      using Base::quadImp;

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

      const typename Base :: Side side_;
      const IntersectionType &intersection_;
      ReferenceGeometry referenceGeometry_;
      mutable CoordinateType dummy_;
    };

    /** \copydoc ElementQuadratureImpl< GridPart, codim, IntegrationTraits > */
    template< class GridPart, int codim, template <class, int> class IntegrationTraits = DefaultQuadratureTraits >
    using ElementIntegrationPointList = ElementQuadratureImpl< GridPart, codim, ElementQuadratureTraits< GridPart, codim, IntegrationTraits>, false >;

    template< class GridPart, int codim, template <class, int> class IntegrationTraits = DefaultQuadratureTraits >
    using ElementQuadrature = ElementQuadratureImpl< GridPart, codim, ElementQuadratureTraits< GridPart, codim, IntegrationTraits>, true >;

    template<class GridPart, class Entity>
    static inline auto elementQuadrature(const GridPart& gridPart, const Entity& entity, unsigned quadOrder)
    {
      using Quadrature = Dune::Fem::ElementQuadrature<GridPart, 0>;
      return Quadrature(entity, quadOrder);
    }

  } // namespace Fem

} // namespace Dune

#include "agglomerationquadrature.hh"

#endif // #ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
