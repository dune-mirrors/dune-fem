#ifndef DUNE_CACHINGPOINTLIST_HH
#define DUNE_CACHINGPOINTLIST_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/quadrature/elementpointlist.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/caching/pointmapper.hh>
#include <dune/fem/quadrature/caching/cacheprovider.hh>

namespace Dune
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
  };



  /** \class CachingPointList
   *  \ingroup Quadrature
   *  \brief integration point list supporting base function caching
   *
   *  A CachingPointList is a conceptual extension to the ElementIntegrationPointList.
   *  It provides an additional mapping from local quadrature point numbers on
   *  a subentity's reference element to global quadrature point numbers on the
   *  codim-0 reference element. Consider, for instance, a quadrature for one
   *  of the faces of a tetrahedron: It provides n local quadrature points, which
   *  can lie on one of the four faces, resulting in 4*n global quadrature points.
   *  
   *  The information from the mapping can be used to cache a base function on
   *  those global quadrature points.
   *
   *  \note If you don't want caching, you can use ElementIntegrationPointList
   *        instead.
   *
   *  For the actual implementation, see
   *  - CachingPointList<GridPartImp,0,IntegrationTraits>
   *  - CachingPointList<GridPartImp,1,IntegrationTraits>
   *
   * \interfaceclass
   */
  template< class GridPartImp, int codim, class IntegrationTraits >
  class CachingPointList;

  
  
  /** \copydoc Dune::CachingPointList */
  template< class GridPartImp, class IntegrationTraits >
  class CachingPointList< GridPartImp, 0, IntegrationTraits >
  : public ElementPointListBase< GridPartImp, 0, IntegrationTraits >,
    public CachingInterface
  {
    typedef CachingPointList< GridPartImp, 0, IntegrationTraits > This;
    typedef ElementPointListBase< GridPartImp, 0, IntegrationTraits > Base;

  public:
    static const int codimension = Base::codimension;

    //! The type of the coordinates in the codim-0 reference element.
    typedef typename Base::CoordinateType CoordinateType;

    //! the type of the quadrature point 
    typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;
    

    // for compatibility
    enum Side { INSIDE, OUTSIDE };
    typedef typename Base::GridPartType::GridType GridType;


  protected:
    using Base::quadImp;

  public:
    using Base::localPoint;

    /** \copydoc Dune::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::ElementIntegrationPointList(const GeometryType &geometry,int order)
     */
    CachingPointList( const GeometryType &geometry, int order )
    : Base( geometry, order )
    {
      CacheProvider< GridType, codimension >::registerQuadrature( quadImp() );
    }

    const QuadraturePointWrapperType operator[] ( const size_t i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }

    /** \copydoc Dune::IntegrationPointList::point */
    const CoordinateType &point ( const size_t i ) const
    {
      return localPoint( i );
    }

    /** \copydoc Dune::CachingInterface::cachingPoint */
    size_t cachingPoint( const size_t i ) const
    {
      return i;
    }
  };
 


  /** \copydoc Dune::CachingPointList */
  template< typename GridPartImp, class IntegrationTraits >
  class CachingPointList< GridPartImp, 1, IntegrationTraits >
  : public ElementPointListBase< GridPartImp, 1, IntegrationTraits >, 
    public CachingInterface
  {
    typedef CachingPointList< GridPartImp, 1, IntegrationTraits > This;
    typedef ElementPointListBase< GridPartImp, 1, IntegrationTraits > Base;

  public:
    //! type of the grid partition
    typedef GridPartImp GridPartType;

    //! side of intersection  
    enum Side { INSIDE, OUTSIDE };

    typedef typename Base::RealType RealType;
    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;

    //! Type of coordinates in codim-0 reference element
    typedef typename Base::CoordinateType CoordinateType;
    
    //! Type of the intersection iterator
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    typedef QuadraturePointWrapper< This > QuadraturePointWrapperType;

    //! type of quadrature used for non-conforming intersections  
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      NonConformingQuadratureType; 


    // for compatibility
    typedef typename GridPartType::GridType GridType;
    typedef TwistUtility< GridType > TwistUtilityType;
    typedef IntersectionIteratorType IntersectionIterator;


  protected:
    typedef typename CachingTraits< RealType, dimension >::MapperType MapperType;
    typedef typename CachingTraits< RealType, dimension >::PointVectorType PointVectorType;

    typedef Dune::CacheProvider< GridType, codimension > CacheProvider;
    typedef Dune::PointProvider< RealType, dimension, codimension> PointProvider;

    using Base::localFaceIndex;
    using Base::quadImp;

  public:
    using Base::elementGeometry;
    using Base::nop;

    /** \brief constructor
     *
     *  \note The CachingPointList requires the grid part to get twist
     *        information for TwistUtility (see also
     *        ElementIntegrationPointList<GridPartImp,1>).
     * 
     *  \param[in]  gridPart      grid partition
     *  \param[in]  intersection  intersection
     *  \param[in]  order         desired order of the quadrature
     *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
     *                            which the ElementQuadrature shall be created
     */
    CachingPointList ( const GridPartType &gridPart,
                       const IntersectionType &intersection,
                       int order, const Side side )
      : Base( getPointList( gridPart, intersection, order, side ) ),
        mapper_( CacheProvider::getMapper( quadImp(), elementGeometry(), localFaceIndex(), twist_ ) ),
        points_( PointProvider::getPoints( quadImp().ipList().id(), elementGeometry() ) )
    {
      //assert( intersection.conforming() );
    }

    const QuadraturePointWrapperType operator[] ( const size_t i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }

    /** \copydoc Dune::IntegrationPointList::point
     */
    const CoordinateType &point ( const size_t i ) const
    {
      return points_[ cachingPoint( i ) ];
    }

    /** \copydoc Dune::CachingInterface::cachingPoint */
    size_t cachingPoint ( const size_t i ) const 
    {
      assert( i < (size_t)nop() );
      return mapper_[ i ];
    }

    // return local caching point 
    // for debugging issues only 
    size_t localCachingPoint ( const size_t i ) const 
    {
      assert( i < (size_t)nop() );

      assert( mapper_[ i ] >= 0 );
      int faceIndex = localFaceIndex();
      int point = mapper_[ i ] - faceIndex * mapper_.size();
      assert( point < nop() );

      return point;
    }

  protected:
    Base getPointList ( const GridPartType &gridPart,
                        const IntersectionType &intersection,
                        const int order, const Side side )
    {
      switch( side )
      {
      case INSIDE:
        twist_ = TwistUtilityType::twistInSelf( gridPart.grid(), intersection );
        return Base( TwistUtilityType::elementGeometry( intersection, true ),
                     intersection.indexInInside(), order );

      case OUTSIDE:
        twist_ = TwistUtilityType::twistInNeighbor( gridPart.grid(), intersection );
        return Base( TwistUtilityType::elementGeometry( intersection, false ),
                     intersection.indexInOutside(), order );

      default:
        DUNE_THROW( InvalidStateException, "ElementIntegrationPointList: side must either be INSIDE or OUTSIDE." );
      }
    }

  private:
    int twist_;
    const MapperType &mapper_;
    const PointVectorType &points_;
  };

}

#endif // #ifndef DUNE_CACHINGPOINTLIST_HH
