#ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
#define DUNE_FEM_ELEMENTQUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/fem/common/utility.hh>

#include "quadrature.hh"
#include "elementpointlist.hh"

namespace Dune
{

  namespace Fem
  {

    template< typename GridPartImp, class IntegrationPointList >
    class Agglomeration;


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
     *  - ElementQuadrature<GridPartImp,0>
     *  - ElementQuadrature<GridPartImp,1>
     */
    template< typename GridPartImp, int codim, template <class, int> class QuadratureTraits = DefaultQuadratureTraits >
    class ElementQuadrature;



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



    /** \copydoc ElementQuadrature */
    template< typename GridPartImp, template< class, int > class QuadratureTraits >
    class ElementQuadrature< GridPartImp, 0, QuadratureTraits >
    : public ElementIntegrationPointList< GridPartImp, 0, ElementQuadratureTraits< GridPartImp, 0, QuadratureTraits > >
    {
      typedef ElementQuadrature< GridPartImp, 0, QuadratureTraits > ThisType;

    public:
      typedef ElementQuadratureTraits< GridPartImp, 0, QuadratureTraits > IntegrationTraits;
      typedef ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits > BaseType;

      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 0 };

      //! dimension of the world
      enum { dimension = GridPartType :: dimension };

      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      //! type for coordinates in the codim-0 reference element
      typedef typename IntegrationTraits :: CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename BaseType :: QuadratureKeyType QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< ThisType > IteratorType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef typename BaseType :: IntegrationPointListType  IntegrationPointListType;

    protected:
      using BaseType::quadImp;

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

    public:
      using BaseType::nop;

      /*! \brief constructor
       *
       *  \param[in]  entity  entity, on whose reference element the quadrature
       *                      lives
       *  \param[in]  quadKey desired minimal order of the quadrature or other means of quadrature identification
       */
      ElementQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey, const bool checkGeomType = true )
      : BaseType( createQuadrature( entity, quadKey, checkGeomType ) )
      {}

      /*! \brief constructor
       *
       *  \param[in]  type    geometry type, on whose reference element the quadrature
       *                      lives
       *  \param[in]  quadKey desired minimal order of the quadrature or other means of quadrature identification
       */
      ElementQuadrature( const GeometryType &type, const QuadratureKeyType& quadKey )
      : BaseType( type, quadKey )
      {
        // when type is none then entity has to be passed in order to create
        // the sub triangulation etc.
        // assert( ! type.isNone() );
      }

      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      ElementQuadrature( const ThisType &org )
      : BaseType( org )
      {}

      QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /** \copydoc Dune::Fem::Quadrature::weight */
      const RealType &weight( size_t i ) const
      {
        return quadImp().weight( i );
      }
    };



    /** \copydoc ElementQuadrature */
    template< class GridPartImp, template< class, int > class QuadratureTraits >
    class ElementQuadrature< GridPartImp, 1, QuadratureTraits >
    : public ElementIntegrationPointList
      < GridPartImp, 1, ElementQuadratureTraits< GridPartImp, 1, QuadratureTraits > >
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the quadrature
      enum { codimension = 1 };

    private:
      typedef ElementQuadratureTraits< GridPartType, codimension, QuadratureTraits > IntegrationTraits;

      typedef ElementQuadrature< GridPartType, codimension, QuadratureTraits > ThisType;
      typedef ElementIntegrationPointList< GridPartType, 1, IntegrationTraits >
        BaseType;

    protected:
      using BaseType :: quadImp;

    public:
      //! dimension of the world
      enum { dimension = GridPartType :: dimension };

      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      //! type of the intersection iterator
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;

      //! type of coordinates in codim-0 reference element
      typedef typename IntegrationTraits :: CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename BaseType :: QuadratureKeyType QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< ThisType > IteratorType;

      //! type of coordinate in codim-1 reference element
      typedef typename IntegrationTraits :: IntegrationPointListType :: CoordinateType
        LocalCoordinateType;

      //! type of quadrature for use on non-conforming intersections
      typedef ThisType NonConformingQuadratureType;

    public:
      using BaseType::nop;

      /*! \brief constructor
       *
       *  \param[in]  gridPart      grid partition (a dummy here)
       *  \param[in]  intersection  intersection
       *  \param[in]  quadKey       quadrature key, i.e. desired order of the quadrature
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for
       *                            which the ElementQuadrature shall be created
       */
      ElementQuadrature ( const GridPartType &gridPart,
                          const IntersectionType &intersection,
                          const QuadratureKeyType& quadKey,
                          typename BaseType :: Side side )
      : BaseType( gridPart, intersection, quadKey, side )
      {}

      /*! \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      ElementQuadrature( const ElementQuadrature &org )
      : BaseType( org )
      {
      }

      QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /*! obtain the weight of the i-th quadrature point
       *
       *  \note The quadrature weights sum up to the volume of the corresponding
       *        reference element.
       *
       *  \param[in]  i  index of the quadrature point
       *
       *  \returns weight of the i-th quadrature point within the quadrature
       */
      const RealType &weight( size_t i ) const
      {
        return quadImp().weight( i );
      }
    };

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
