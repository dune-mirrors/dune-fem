#ifndef DUNE_FEM_CACHINGQUADRATURE_HH
#define DUNE_FEM_CACHINGQUADRATURE_HH

//- Dune includes
#include <dune/common/math.hh>

//- Local includes
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/caching/pointmapper.hh>
#include <dune/fem/quadrature/caching/cacheprovider.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class CachingInterface
     *  \brief interface a cachable quadrature has to implement
     */
    class CachingInterface
    {
    protected:
      // do not create instances of this class
      CachingInterface ()
      {
      }

    public:
      /** \brief returns true if cachingPoint is not the identity mapping */
      static constexpr bool twisted () { return false; }

      /** \brief returns the twistId, i.e. [0,...,7] */
      inline int twistId () const { return 0; }

      /** \brief map quadrature points to caching points
       *
       *  For codim-1 entites, the mapping consists of two stages:
       *  - Consider the twist to get the quadrature point number on the face of
       *    the (codim-0) reference element,
       *  - Map the twisted quadrature point number to the caching point number.
       *
       *  \param[in]  quadraturePoint  number of quadrature point to map to a
       *                               caching point
       */
      inline size_t cachingPoint( const size_t quadraturePoint ) const
      {
        DUNE_THROW( NotImplemented,
                    "CachingInterface :: cachingPoint must be overloaded!" );
      }

      /** \brief map quadrature points to interpolation points
       *
       *  \param[in]  quadraturePoint  number of quadrature point to map to an
       *                               interpolation point
       */
      inline size_t interpolationPoint( const size_t quadraturePoint ) const
      {
        DUNE_THROW( NotImplemented,
                    "CachingInterface :: interpolationPoint must be overloaded!" );
      }

      /** \brief check if quadrature is interpolation quadrature
       *
       *  \param[in]  numShapeFunctions  number of shapeFunctions that has to
       *              match number of quadrature points or number of
       *              internal interpolation  points
       */
      inline bool isInterpolationQuadrature( const size_t numShapeFunctions ) const
      {
        DUNE_THROW( NotImplemented,
                    "CachingInterface :: isInterpolationQuadrature must be overloaded!" );
      }
    };



    /** \class CachingQuadratureImpl
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
     * \interfaceclass
     */
    template< class GridPartImp, int codim, class IntegrationTraits, bool isQuadrature >
    class CachingQuadratureImpl;


    /** \copydoc CachingQuadratureImpl  */
    template< class GridPartImp, class IntegrationTraits, bool isQuadrature >
    class CachingQuadratureImpl< GridPartImp, 0, IntegrationTraits, isQuadrature >
    : public ElementPointListBase< GridPartImp, 0, IntegrationTraits >,
      public CachingInterface
    {
      typedef CachingQuadratureImpl< GridPartImp, 0, IntegrationTraits, isQuadrature > This;
      typedef ElementPointListBase< GridPartImp, 0, IntegrationTraits > Base;

    public:
      typedef typename Base::GridPartType GridPartType;
      typedef typename Base::EntityType   EntityType;
      static const int codimension = Base::codimension;

      //! The type of the coordinates in the codim-0 reference element.
      typedef typename Base::CoordinateType CoordinateType;

      //! the type of the quadrature point
      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< This > IteratorType;

      //! id of point set, positive if interpolation point set, otherwise negative
      static const int pointSetId = SelectQuadraturePointSetId<
          typename IntegrationTraits::IntegrationPointListType::Traits > :: value;

    protected:
      using Base::quadImp;

    public:
      using CachingInterface::twisted;
      using CachingInterface::twistId;
      using Base::localPoint;
      using Base::nop;

      /** \copydoc Dune::Fem::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::ElementIntegrationPointList(const GeometryType &geometry, const QuadratureKeyType& quadKey)
       */
      template <class QuadratureKeyType>
      CachingQuadratureImpl( const EntityType& entity, const QuadratureKeyType& quadKey )
      : Base( entity.type(), quadKey )
      {
        CacheProvider< GridPartType, codimension >::registerQuadrature( quadImp() );
      }

      /** \copydoc Dune::Fem::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::ElementIntegrationPointList(const GeometryType &geometry, const QuadratureKeyType& quadKey)
       */
      template <class QuadratureKeyType>
      CachingQuadratureImpl( const GeometryType &geometry, const QuadratureKeyType& quadKey )
      : Base( geometry, quadKey )
      {
        CacheProvider< GridPartType, codimension >::registerQuadrature( quadImp() );
      }

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
      auto weight ( std::size_t i ) const
      {
        return quadImp().weight( i );
      }

      /** \copydoc Dune::Fem::CachingInterface::cachingPoint */
      inline size_t cachingPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      /** \copydoc Dune::Fem::CachingInterface::interpolationPoint */
      inline size_t interpolationPoint( const size_t quadraturePoint ) const
      {
        return quadraturePoint;
      }

      /** \copydoc Dune::Fem::CachingInterface::isInterpolationQuadrature */
      inline bool isInterpolationQuadrature( const size_t numShapeFunctions ) const
      {
        // if pointSetId is not negative then we have an interpolation
        // quadrature if the number of point are equal to number of shape functions
        return (pointSetId >= 0) ? (nop() == numShapeFunctions) : false;
      }
    };



    /** \copydoc CachingQuadratureImpl */
    template< typename GridPartImp, class IntegrationTraits, bool isQuadrature >
    class CachingQuadratureImpl< GridPartImp, 1, IntegrationTraits, isQuadrature >
    : public ElementPointListBase< GridPartImp, 1, IntegrationTraits >,
      public CachingInterface
    {
      typedef CachingQuadratureImpl< GridPartImp, 1, IntegrationTraits, isQuadrature > This;
      typedef ElementPointListBase< GridPartImp, 1, IntegrationTraits > Base;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      typedef typename Base::RealType RealType;
      static const int dimension = Base::dimension;
      static const int codimension = Base::codimension;

      //! Type of coordinates in codim-0 reference element
      typedef typename Base::CoordinateType CoordinateType;

      //! Type of the intersection iterator
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< This > IteratorType;

      //! type of quadrature used for non-conforming intersections
      typedef ElementQuadratureImpl< GridPartType, codimension, IntegrationTraits, isQuadrature >
        NonConformingQuadratureType;


      // for compatibility
      typedef typename GridPartType::TwistUtilityType  TwistUtilityType;
      typedef IntersectionIteratorType IntersectionIterator;

    private:
      static const int quadPointSetId =
        SelectQuadraturePointSetId< typename IntegrationTraits::IntegrationPointListType::Traits > :: value;

    public:
      // Note: we also exclude GaussLegendre(0) here, because on faces it is not
      //       an interpolation rule
      static const int pointSetId = (quadPointSetId > 0) ? quadPointSetId :
                  SelectQuadraturePointSetId< void > :: value; // default value

    protected:
      typedef typename CachingTraits< RealType, dimension >::MapperPairType  MapperPairType;
      typedef typename CachingTraits< RealType, dimension >::PointVectorType PointVectorType;

      typedef CacheProvider< GridPartType, codimension >            CacheProviderType;
      typedef PointProvider< RealType, dimension, codimension>  PointProviderType;

      using Base::quadImp;

    public:
      using Base::localFaceIndex;
      using Base::elementGeometry;
      using Base::nop;

      /** \brief constructor
       *
       *  \note The CachingQuadratureImpl requires the grid part to get twist
       *        information for TwistUtility (see also
       *        ElementIntegrationPointList<GridPartImp,1>).
       *
       *  \param[in]  gridPart      grid partition
       *  \param[in]  intersection  intersection
       *  \param[in]  quadKey       desired order of the quadrature or other means of quadrature identification
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for
       *                            which the ElementQuadrature shall be created
       */
      template <class QuadratureKeyType>
      CachingQuadratureImpl ( const GridPartType &gridPart,
                              const IntersectionType &intersection,
                              const QuadratureKeyType& quadKey, const typename Base :: Side side )
        : Base( getPointList( intersection, quadKey, side ) ),
          side_(side),
          twist_( getTwist( gridPart, intersection, side ) ),
          mapper_( CacheProviderType::getMapper( quadImp(), elementGeometry(), localFaceIndex(), twist_) ),
          points_( PointProviderType::getPoints( quadImp().ipList().id(), elementGeometry() ) ),
          intersection_(intersection)
      {
      }

      const QuadraturePointWrapperType operator[] ( const size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, nop() ); }

      typename Base :: Side side() const { return side_; }
      bool isInside() const { return side_ == Base::INSIDE; }

      /** \copydoc Dune::Fem::IntegrationPointList::point
       */
      const CoordinateType &point ( const size_t i ) const
      {
        return points_[ cachingPoint( i ) ];
      }

      /** \copydoc Dune::Fem::IntegrationPointList::weight */
      auto weight ( std::size_t i ) const
      {
        return quadImp().weight( i );
      }

      const IntersectionType &intersection() const
      {
        return intersection_;
      }

      /** \copydoc Dune::Fem::CachingInterface::twisted */
      static constexpr bool twisted() { return true; }

      /** \copydoc Dune::Fem::CachingInterface::twistId */
      inline int twistId () const { return twist_ + 4; }

      /** \copydoc Dune::Fem::CachingInterface::cachingPoint */
      inline size_t cachingPoint( const size_t quadraturePoint ) const
      {
        assert( quadraturePoint < (size_t)nop() );
        return mapper_.first[ quadraturePoint ];
      }

      /** \copydoc Dune::Fem::CachingInterface::interpolationPoint */
      inline size_t interpolationPoint( const size_t quadraturePoint ) const
      {
        assert( quadraturePoint < mapper_.second.size() );
        return mapper_.second[ quadraturePoint ];
      }

      /** \copydoc Dune::Fem::CachingInterface::isInterpolationQuadrature */
      inline bool isInterpolationQuadrature( const size_t numShapeFunctions ) const
      {
        // if pointSetId is not negative then we have an interpolation
        // quadrature if the number of point are equal to number of shape functions
        return (pointSetId < 0) ? false :
          quadImp().ipList().isFaceInterpolationQuadrature( numShapeFunctions );
      }

      // return local caching point
      // for debugging issues only
      size_t localCachingPoint ( const size_t i ) const
      {
        const auto& mapper = mapper_.first;

        assert( i < (size_t)nop() );

        assert( mapper[ i ] >= 0 );
        int faceIndex = localFaceIndex();
        unsigned int point = mapper[ i ] - faceIndex * mapper.size();
        assert( point < nop() );

        return point;
      }

    protected:
      template <class QuadratureKeyType>
      Base getPointList ( const IntersectionType &intersection,
                          const QuadratureKeyType& key,
                          const typename Base :: Side side )
      {
        switch( side )
        {
          case Base :: INSIDE:
            return Base( TwistUtilityType::elementGeometry( intersection, true ),
                         intersection.indexInInside(), key );

          case Base :: OUTSIDE:
            return Base( TwistUtilityType::elementGeometry( intersection, false ),
                         intersection.indexInOutside(), key );

          default:
            DUNE_THROW( InvalidStateException, "ElementIntegrationPointList: side must either be INSIDE or OUTSIDE." );
        }
      }

      int getTwist ( const GridPartType &gridPart,
                     const IntersectionType &intersection,
                     const typename Base :: Side side )
      {
        switch( side )
        {
          case Base :: INSIDE:
            return TwistUtilityType::twistInSelf( gridPart.grid(), intersection );

          case Base :: OUTSIDE:
            return TwistUtilityType::twistInNeighbor( gridPart.grid(), intersection );

          default:
            DUNE_THROW( InvalidStateException, "ElementIntegrationPointList: side must either be INSIDE or OUTSIDE." );
        }
      }

    private:
      const typename Base :: Side side_;
      const int twist_;
      const MapperPairType &mapper_;
      const PointVectorType &points_;
      const IntersectionType &intersection_;
    };


    /** \copydoc CachingQuadratureImpl */
    template< typename GridPart, int codim, class IntegrationTraits >
    using CachingPointList = CachingQuadratureImpl< GridPart, codim, IntegrationTraits, false>;

    /** \copydoc CachingQuadratureImpl */
    template< typename GridPart, int codim, template <class, int> class QuadratureTraits = DefaultQuadratureTraits >
    using CachingQuadrature = CachingQuadratureImpl< GridPart, codim, ElementQuadratureTraits< GridPart, codim, QuadratureTraits >, true>;

    template<class GridPart, class Entity>
    static inline auto cachingQuadrature(const GridPart& gridPart, const Entity& entity, unsigned quadOrder)
    {
      using Quadrature = Dune::Fem::CachingQuadrature<GridPart, 0>;
      return Quadrature(entity, quadOrder);
    }

  } //namespace Fem

} //namespace Dune

#endif // #ifndef DUNE_FEM_CACHINGQUADRATURE_HH
