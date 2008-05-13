#ifndef DUNE_CACHEPOINTLIST_HH
#define DUNE_CACHEPOINTLIST_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include <dune/fem/quadrature/elementquadrature.hh>
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
  class CachingPointList
  {
    typedef CompileTimeChecker< false > Only_specialisations_for_codim_0_and_1_so_far;
  };

  
  
  /** \copydoc Dune::CachingPointList */
  template< class GridPartImp, class IntegrationTraits >
  class CachingPointList< GridPartImp, 0, IntegrationTraits >
  : public ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits >,
    public CachingInterface
  {
  public:
    //! type of grid partition
    typedef GridPartImp GridPartType;

    //! codimension of element quadrature
    enum { codimension = 0 };
    
  private:
    typedef CachingPointList< GridPartType, codimension, IntegrationTraits > ThisType;
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      BaseType;

  protected:
    using BaseType :: quadImp;

  public:
    // type of grid 
    typedef typename GridPartImp :: GridType GridType;

    //! Dimension of the world.
    enum { dimension = BaseType::dimension };

    //! Just another name for double...
    typedef typename BaseType::RealType RealType;
    
    //! The type of the coordinates in the codim-0 reference element.
    typedef typename BaseType::CoordinateType CoordinateType;

    //! The type of the codim-0 entity.
    typedef typename BaseType::Entity Entity;

    //! the type of the quadrature point 
    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
    
  public:
    /** \copydoc Dune::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::ElementIntegrationPointList(const GeometryType &geometry,int order)
     */
    inline CachingPointList( const GeometryType &geometry, int order )
    : BaseType( geometry, order )
    {
      CacheProvider< GridType, codimension > :: registerQuadrature( quadImp() );
    }

    /** \brief copy constructor
     *
     *  \param[in]  org  element quadrature to copy
     */
    inline CachingPointList( const ThisType& org )
    : BaseType( org )
    {
    }

    inline const QuadraturePointWrapperType operator[] ( const unsigned int i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }

    /** \copydoc Dune::CachingInterface::cachingPoint */
    inline size_t cachingPoint( const size_t quadraturePoint ) const
    {
      return quadraturePoint;
    }
  };
 


  /** \copydoc Dune::CachingPointList */
  template< typename GridPartImp, class IntegrationTraits >
  class CachingPointList< GridPartImp, 1, IntegrationTraits >
  : public ElementIntegrationPointList< GridPartImp, 1, IntegrationTraits >, 
    public CachingInterface
  {
  public:
    //! type of grid partition
    typedef GridPartImp GridPartType;

    //! codimension of the element quadrature
    enum { codimension = 1 };
    
  private:
    typedef CachingPointList< GridPartType, codimension, IntegrationTraits > ThisType;
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      BaseType;

 protected:
    using BaseType :: faceNumber;
    using BaseType :: quadImp;

  public:
    using BaseType :: elementGeometry;
    using BaseType :: nop;
    
  public:
    //! type of the grid
    typedef typename GridPartType :: GridType GridType;

    //! Dimeinsion of the world
    enum { dimension = BaseType::dimension };
    
    //! A double... or whatever your grid wants
    typedef typename BaseType::RealType RealType;
    
    //! The coordinates of the quadrature points in the codim-0 reference
    //! element
    typedef typename BaseType::CoordinateType CoordinateType;

    //! Type of the intersection iterator
    typedef typename BaseType::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection IntersectionType;

    //! type of quadrature used for non-conforming intersections  
    typedef BaseType NonConformingQuadratureType; 

    //! type of twist utility 
    typedef TwistUtility< GridType > TwistUtilityType;

    //! the type of the quadrature point 
    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
    
  protected:
    typedef typename CachingTraits< RealType, dimension > :: MapperType MapperType;

  protected:
    const MapperType &mapper_;

  public:
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
                       int order,
                       typename BaseType :: Side side)
      : BaseType( gridPart, intersection, order, side ),
        mapper_( CacheProvider< GridType, codimension > :: getMapper
          ( quadImp(), elementGeometry(), faceNumber(),
            (side == BaseType :: INSIDE)
              ? TwistUtilityType :: twistInSelf( gridPart.grid(), intersection )
              : TwistUtilityType :: twistInNeighbor( gridPart.grid(), intersection )
          ) )
    {
      // make sure CachingPointList is only created for conforming intersections
      // assert( TwistUtilityType :: conforming( gridPart.grid(), intersection ) );
    }

    /** \brief copy constructor
     *
     *  \param[in]  org  element quadrature to copy
     */
    CachingPointList( const ThisType& org )
    : BaseType( org ),
      mapper_( org.mapper_ )
    {
    }

    inline const QuadraturePointWrapperType operator[] ( const unsigned int i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }

    /** \copydoc Dune::CachingInterface::cachingPoint */
    size_t cachingPoint( const size_t quadraturePoint ) const 
    {
      assert( quadraturePoint < (size_t)nop() );
      return mapper_[ quadraturePoint ];
    }

    // return local caching point 
    // for debugging issues only 
    size_t localCachingPoint(size_t quadraturePoint) const 
    {
      // this makes no sense for usigned ints ;)
      assert(quadraturePoint >= 0);
      assert(quadraturePoint < (size_t)this->nop());

      int faceIndex = this->faceNumber();
      int point = mapper_[quadraturePoint] - faceIndex*mapper_.size();
      assert( mapper_[quadraturePoint] >= 0 );

      assert( point < this->nop() );
      return point;
    }
  };
}

#endif
