#ifndef DUNE_FEM_CACHINGQUADRATURE_HH
#define DUNE_FEM_CACHINGQUADRATURE_HH

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

#include "cachingpointlist.hh"

namespace Dune
{
  namespace Fem
  {

    /** \class CachingQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature class supporting base function caching
     *
     *  A CachingQuadrature is a conceptual extension to the ElementQuadrature.
     *  It provides an additional mapping from local quadrature point numbers on
     *  a subentity's reference element to global quadrature point numbers on the
     *  codim-0 reference element. Consider, for instance, a quadrature for one
     *  of the faces of a tetrahedron: It provides n local quadrature points, which
     *  can lie on one of the four faces, resulting in 4*n global quadrature points.
     *
     *  The information from the mapping can be used to cache a base function on
     *  those global quadrature points.
     *
     *  \note If you don't want caching, you can use ElementQuadrature instead.
     *
     *  \note QuadratureTraits specifies the quadrature points used to build the
     *        quadrature.
     *
     *  For the actual implementations, see
     *  - CachingQuadrature<GridPartImp,0>
     *  - CachingQuadrature<GridPartImp,1>
     */
    template< typename GridPartImp, int codim, template <class, int> class QuadratureTraits = DefaultQuadratureTraits >
    class CachingQuadrature;



    /** \copydoc CachingQuadrature */
    template< typename GridPart, template <class, int> class QuadratureTraits >
    class CachingQuadrature< GridPart, 0, QuadratureTraits >
    : public CachingPointList< GridPart, 0, ElementQuadratureTraits< GridPart, 0, QuadratureTraits > >
    {
    public:
      //! type of grid partition
      typedef GridPart GridPartType;

      //! codimension of the element quadrature
      static constexpr auto codimension = 0;

    private:
      typedef ElementQuadratureTraits< GridPartType, codimension, QuadratureTraits > IntegrationTraits;

      typedef CachingQuadrature< GridPartType, codimension, QuadratureTraits > ThisType;
      typedef CachingPointList< GridPartType, codimension, IntegrationTraits > BaseType;

    public:
      //! dimension of the world
      static constexpr auto dimension = BaseType::dimension;

      //! just another name for double
      typedef typename BaseType :: RealType RealType;
      //! type of the coordinates in the codim-0 reference element
      typedef typename BaseType :: CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename BaseType::QuadratureKeyType QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< ThisType > IteratorType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    protected:
      using BaseType :: quadImp;

    public:
      using BaseType::nop;

      /** \brief constructor
       *
       *  \param[in]  entity  entity, on whose reference element the quadrature
       *                      lives
       *  \param[in]  quadKey desired minimal order of the quadrature or other means of quadrature identification
       */
      CachingQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey )
      : BaseType( entity.type(), quadKey )
      {}

      /** \brief constructor
       *
       *  \param[in]  type    geometry type, on whose reference element the quadrature
       *                      lives
       *  \param[in]  quadKey desired minimal order of the quadrature or other means of quadrature identification
       */
      CachingQuadrature( const GeometryType &type, const QuadratureKeyType& quadKey )
      : BaseType( type, quadKey )
      {}

      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      CachingQuadrature( const ThisType &org )
      : BaseType( org )
      {}

      QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /** \copydoc Dune::Fem::ElementQuadrature<GridPartImp,0>::weight */
      const RealType &weight ( std::size_t i ) const
      {
        return quadImp().weight( i );
      }
    };



    /** \copydoc CachingQuadrature */
    template< typename GridPartImp, template <class, int> class QuadratureTraits  >
    class CachingQuadrature< GridPartImp, 1, QuadratureTraits >
    : public CachingPointList
      < GridPartImp, 1, ElementQuadratureTraits< GridPartImp, 1, QuadratureTraits > >
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      static constexpr auto codimension = 1;

    private:
      typedef ElementQuadratureTraits< GridPartType, codimension, QuadratureTraits > IntegrationTraits;

      typedef CachingQuadrature< GridPartType, codimension, QuadratureTraits > ThisType;
      typedef CachingPointList< GridPartType, codimension, IntegrationTraits > BaseType;

    protected:
      using BaseType :: quadImp;

    public:
      //! dimeinsion of the world
      static constexpr auto dimension = BaseType::dimension;

      //! just another name for double
      typedef typename BaseType::RealType RealType;

      //! the coordinates of the quadrature points in the codim-0 reference element
      typedef typename BaseType::CoordinateType CoordinateType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef typename BaseType::QuadratureKeyType QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< ThisType > IteratorType;

      //! type of the intersection iterator
      typedef typename BaseType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;

      //! type of quadrature used for non-conforming intersections
      typedef ElementQuadrature< GridPartImp, codimension > NonConformingQuadratureType;

      using BaseType::nop;

      /** \brief constructor
       *
       *  \note The CachingQuadrature requires the grid part to get twist
       *        information for TwistUtility (see also
       *        ElementQuadrature<GridPartImp,1>).
       *
       *  \param[in]  gridPart      grid partition
       *  \param[in]  intersection  intersection
       *  \param[in]  order         desired order of the quadrature
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for
       *                            which the ElementQuadrature shall be created
       */
      CachingQuadrature( const GridPartType &gridPart, const IntersectionType &intersection,
                         const QuadratureKeyType& quadKey, typename BaseType::Side side )
      : BaseType( gridPart, intersection, quadKey, side )
      {}

      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      CachingQuadrature( const ThisType& org )
      : BaseType( org )
      {}

      QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      /** \copydoc Dune::Fem::ElementQuadrature<GridPartImp,1>::weight */
      const RealType &weight( std::size_t i ) const
      {
        return quadImp().weight(i);
      }
    };

    template<class GridPart, class Entity>
    static inline auto cachingQuadrature(const GridPart& gridPart, const Entity& entity, unsigned quadOrder)
    {
      using Quadrature = Dune::Fem::CachingQuadrature<GridPart, 0>;
      return Quadrature(entity, quadOrder);
    }
  } //namespace Fem

} //namespace Dune

#endif // #ifndef DUNE_FEM_CACHINGQUADRATURE_HH
