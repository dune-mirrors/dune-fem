#ifndef DUNE_ELEMENTPOINTLIST_HH
#define DUNE_ELEMENTPOINTLIST_HH

#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
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
  class ElementIntegrationPointList
  {
    typedef CompileTimeChecker< false >
      __Only_implementations_for_codim_0_and_1_exist__; 
  };



  /** \copydoc Dune::ElementIntegrationPointList */
  template< class GridPartImp, class IntegrationTraits >
  class ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits >
  {
  public:
    //! type of the grid partition
    typedef GridPartImp GridPartType;

    //! codimension of the integration point list
    enum { codimension = 0 };

  private:
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      ThisType;

  public:
    //! type of grid
    typedef typename GridPartType :: GridType GridType;
    
    //! dimension of the world
    enum { dimension = GridType :: dimension };

     //! coordinate type 
    typedef typename GridType :: ctype RealType;
   
    //! side of intersection 
    enum Side { INSIDE, OUTSIDE };

    //! type of the integration point list 
    typedef typename IntegrationTraits ::  IntegrationPointListType
      IntegrationPointListType;

    //! type for coordinates in the codim-0 reference element
    typedef typename IntegrationTraits :: CoordinateType CoordinateType;

    //! type of the codim-0 entity
    typedef typename GridType :: template Codim< 0 > :: Entity Entity;
    
    //! the type of the quadrature point 
    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
    
  protected:
    IntegrationPointListType quad_;
   
  public:
    /** \brief constructor
     *  
     *  \param[in]  geometry  geometry type, the quadrature lives on
     *  \param[in]  order     desired minimal order of the quadrature
     */
    ElementIntegrationPointList( const GeometryType &geometry, int order )
    : quad_( geometry, order )
    {
    }
    
    /** \brief copy constructor
     *
     *  \param[in]  org  element integration point list to copy
     */
    ElementIntegrationPointList( const ElementIntegrationPointList &org )
    : quad_( org.quad_ )
    {
    }

    inline const QuadraturePointWrapperType operator[] ( unsigned int i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }
   
    /** \copydoc Dune::IntegrationPointList::nop */
    inline int nop () const
    {
      return quad_.nop();
    }

    /** \copydoc Dune::IntegrationPointList::point */
    inline const CoordinateType &point ( size_t i ) const
    {
      return quad_.point(i);
    }
    
    /** \brief obtain local coordinates of i-th integration point
     *
     *  This method returns a reference to the local coordinates of the i-th
     *  integration point for 0 <= i < nop(). Here, local coordinates means
     *  coordinates with respect to the reference element of the subentity.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns reference to i-th integration point
     */
    const CoordinateType &localPoint( size_t i ) const
    {
      return quad_.point(i);
    }

    /** \copydoc Dune::IntegrationPointList::id
     */
    size_t id () const
    {
      return quad_.id();
    }

    /** \copydoc Dune::IntegrationPointList::order
     */
    int order () const
    {
      return quad_.order();
    }

    /** \copydoc Dune::IntegrationPointList::geometry
     */
    GeometryType geometry () const
    {
      return quad_.geometry();
    }
    
    /** \brief obtain GeometryType of the corresponding codim-0 the integration
     *         point list belongs to
     *
     *  An element integration point list can return the coordinates of integration
     *  points with resepct to the codim-0 reference element and the reference
     *  element corresponding to the subentity the quadrature actually lives on.
     *  This method returns the geometry of the codim-0 entity.
     *
     *  \note Calling this method yields a virtual function call, so do not
     *        call this method unnecessarily.
     *
     *  \returns GeometryType for this integration point list
     */
    GeometryType elementGeometry () const
    {
      return quad_.geometry();
    }

  protected:
    /** \brief obtain the actual implementation of the quadrature
     *
     *  \note This method may only be used in derived classes.
     *
     *  \returns a reference to the actual implementation of the quadrature
     */
    const IntegrationPointListType &quadImp () const
    {
      return quad_;
    }
  };


  /** \copydoc Dune::ElementIntegrationPointList */
  template< class GridPartImp, class IntegrationTraits >
  class ElementIntegrationPointList< GridPartImp, 1, IntegrationTraits >
  {
  public:
    //! type of the grid partition
    typedef GridPartImp GridPartType;

    //! codimension of the element integration point list
    enum { codimension = 1 };

  private:
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      ThisType;

  public:
    //! type of the grid 
    typedef typename GridPartType :: GridType GridType;

    //! dimension of the world
    enum { dimension = GridType::dimension };
    
    //! side of intersection  
    enum Side { INSIDE, OUTSIDE };

    //! coordinate type 
    typedef typename GridType :: ctype RealType;

    //! type of the integration point list 
    typedef typename IntegrationTraits :: IntegrationPointListType
      IntegrationPointListType;

    //! Type of coordinates in codim-0 reference element
    typedef typename IntegrationTraits :: CoordinateType CoordinateType;
    
    //! Type of coordinate in codim-1 reference element
    typedef typename IntegrationPointListType :: CoordinateType LocalCoordinateType;

    //! Type of the intersection iterator
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    // For compatibility
    typedef IntersectionIteratorType IntersectionIterator;

    //! type quadrature for use on non-conforming intersections 
    typedef ThisType NonConformingQuadratureType;
   
    //! the type of the quadrature point 
    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;

    //! type of twist utility 
    typedef TwistUtility< GridType > TwistUtilityType;

  private:
    typedef typename IntersectionIteratorType :: LocalGeometry ReferenceGeometry;

  protected:
    const ReferenceGeometry &referenceGeometry_;
    const GeometryType elementGeometry_;
    const IntegrationPointListType quad_;
    const int faceNumber_;

    mutable CoordinateType dummy_;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridPart      grid partition (a dummy here)
     *  \param[in]  intersection  intersection
     *  \param[in]  order         desired order of the quadrature
     *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
     *                            which the ElementQuadrature shall be created
     *
     *  \note This code assumes that the codim-0 entity is either a simplex or
     *        a cube (otherwise elementGeometry() returns a wrong geometry).
     */
    ElementIntegrationPointList ( const GridPartType &gridPart, 
                                  const IntersectionType &intersection, 
                                  int order,
                                  Side side )
    : referenceGeometry_( side == INSIDE ? intersection.intersectionSelfLocal() 
                                         : intersection.intersectionNeighborLocal() ),
      elementGeometry_( TwistUtilityType::elementGeometry(intersection, side == INSIDE ) ), 
      quad_( referenceGeometry_.type() , order ),
      faceNumber_( side == INSIDE ? intersection.numberInSelf()
                                  : intersection.numberInNeighbor() ),
      dummy_( 0. )
    {
    }
    
    /** \brief copy constructor
     *
     *  \param[in]  org  element quadrature to copy
     */
    ElementIntegrationPointList ( const ElementIntegrationPointList &org )
    : referenceGeometry_( org.referenceGeometry_ ),
      elementGeometry_( org.elementGeometry_ ),
      quad_( org.quad_ ),
      faceNumber_( org.faceNumber_ ),
      dummy_( org.dummy_ )
    {
    }
    
    inline const QuadraturePointWrapperType operator[] ( size_t i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }
   
    /** \copydoc Dune::IntegrationPointList::nop
     */
    int nop () const
    {
      return quad_.nop();
    }

    /** \copydoc Dune::IntegrationPointList::point
     */
    const CoordinateType &point ( size_t i ) const
    {
      dummy_ = referenceGeometry_.global( quad_.point( i ) );
      return dummy_;
    }

    /** \copydoc Dune::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::localPoint(size_t i) const */
    const LocalCoordinateType &localPoint ( size_t i ) const
    {
      return quad_.point( i );
    }

    /** \copydoc Dune::IntegrationPointList::id
     */
    size_t id () const
    {
      return quad_.id();
    }

    /** \copydoc Dune::IntegrationPointList::order
     */
    int order () const
    {
      return quad_.order();
    }

    /** \copydoc Dune::IntegrationPointList::geometry
     */
    GeometryType geometry () const
    {
      return quad_.geo();
    }

    /** \copydoc Dune::ElementIntegrationPointList<GridPartImp,0,IntegrationTraits>::elementGeometry() const */
    const GeometryType& elementGeometry () const
    {
      return elementGeometry_;
    }

  protected:
    // return local face number 
    int faceNumber() const
    {
      return faceNumber_;
    }

    /** \brief obtain the actual implementation of the quadrature
     *
     *  \note This method may only be used in derived classes.
     *
     *  \returns a reference to the actual implementation of the quadrature
     */
    const IntegrationPointListType &quadImp() const
    {
      return quad_;
    }
  };

} // end namespace Dune
#endif
