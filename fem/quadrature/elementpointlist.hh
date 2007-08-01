#ifndef DUNE_ELEMENTPOINTLIST_HH
#define DUNE_ELEMENTPOINTLIST_HH

#include "quadrature.hh"

namespace Dune
{
  
  /*! \class ElementIntegrationPointList
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
   *  reference element. Moreover, the transformation is done only once (when
   *  the quadrature is created).
   *  
   *  To achieve this goal, an ElementIntegrationPointList depends stronger on
   *  the context in which it is used. For example, for each face within a
   *  tetrahedron (though they are all the same) we need a different
   *  ElementIntegrationPointList, since the coordinates of the quadrature points
   *  with respect to the codim-0 entity differ for each face.
   *
   *  \note Actually, codim-1 element integration point lists depend on the
   *        intersection.
   */
  template< class GridPartImp, int codim, class IntegrationTraits >
  class ElementIntegrationPointList
  {
    typedef CompileTimeChecker< false >
      __Only_implementations_for_codim_0_and_1_exist__; 
  };



  // \brief Element quadrature on codim-0 entities.
  // For codim-0 element quadratures, there is no additional information
  // from the context needes, in consequence, the quadrature behaves like
  // a generic quadrature class, independent from the situation in the grid.
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
    
  protected:
    IntegrationPointListType quad_;
   
  public:
    /*! \brief constructor
     *  
     *  \param[in]  geometry  geometry type, the quadrature lives on
     *  \param[in]  otder     desired minimal order of the quadrature
     */
    ElementIntegrationPointList( const GeometryType &geometry, int order )
    : quad_( geometry, order )
    {
    }
    
    /*! \brief copy constructor
     *
     *  \param[in]  org  element integration point list to copy
     */
    ElementIntegrationPointList( const ElementIntegrationPointList &org )
    : quad_( org.quad_ )
    {
    }
   
    //! \copydoc Dune::IntegrationPointList::nop
    int nop () const
    {
      return quad_.nop();
    }

    //! \copydoc Dune::IntegrationPointList::point
    const CoordinateType &point ( size_t i ) const
    {
      return quad_.point(i);
    }
    
    /*! \brief obtain local coordinates of i-th integration point
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

    //! \copydoc Dune::IntegrationPointList::id
    size_t id () const
    {
      return quad_.id();
    }

    //! \copydoc Dune::IntegrationPointList::order
    int order () const
    {
      return quad_.order();
    }

    //! \copydoc Dune::IntegrationPointList::geometry
    GeometryType geometry () const
    {
      return quad_.geometry();
    }
    
    /*! \brief obtain GeometryType of the corresponding codim-0 the integration
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
    /*! \brief obtain the actual implementation of the quadrature
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


  //! \brief Element quadrature on codim-1 entities.
  //! For codimension 1, the quadrature needs information about the 
  //! intersection. Plus, the user must decide if the quadrature shall live
  //! on the reference element of the outside or inside element of the 
  //! intersection.
  template <class GridPartImp, class IntegrationTraits>
  class ElementIntegrationPointList<GridPartImp,1,IntegrationTraits>
  {
    typedef GridPartImp GridPartType;
    typedef typename GridPartType :: GridType GridType;

    //! type of this class 
    typedef ElementIntegrationPointList<GridPartImp,1,IntegrationTraits> ThisType;
  public:
    //! Dimension of the world
    enum { dimension = GridType::dimension };
    //! The codimension is one by definition
    enum { codimension = 1 };
    
    //! side of intersection  
    enum Side { INSIDE, OUTSIDE };

    //! coordinate type 
    typedef typename GridType::ctype RealType;

    //! type of the integration point list 
    typedef typename IntegrationTraits :: 
          IntegrationPointListType   IntegrationPointListType;

    //! Type of coordinates in codim-0 reference element
    typedef typename IntegrationTraits :: CoordinateType CoordinateType;
    
    //! Type of coordinate in codim-1 reference element
    typedef typename IntegrationPointListType::CoordinateType LocalCoordinateType;

    //! Type of the intersection iterator
    typedef typename GridPartImp::IntersectionIteratorType IntersectionIterator;

    //! specify quadrature for use on conforming and non-conforming
    //! intersections 
    typedef ThisType NonConformingQuadratureType;
    
  public:
    //! Constructor
    //! \param gridPart s dummy parameter here 
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param side Is either INSIDE or OUTSIDE
    ElementIntegrationPointList(const GridPartType & gridPart, 
                                const IntersectionIterator& it, 
                                int order, Side side) :
      quad_(it.intersectionGlobal().type(), order),
      referenceGeometry_(side == INSIDE ?
                         it.intersectionSelfLocal() : 
                         it.intersectionNeighborLocal()),
      elementGeometry_(referenceGeometry_.type().basicType() ,dimension),
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      dummy_(0.)
    {
    }

    //! Constructor
    //! \param it Intersection iterator
    //! \param order Desired order of the quadrature
    //! \param side Is either INSIDE or OUTSIDE
    ElementIntegrationPointList(const IntersectionIterator& it, int order, Side side) :
      quad_(it.intersectionGlobal().type(), order),
      referenceGeometry_(side == INSIDE ?
                         it.intersectionSelfLocal() : 
                         it.intersectionNeighborLocal()),
      elementGeometry_(referenceGeometry_.type().basicType() ,dimension),
      faceNumber_(side == INSIDE ?
                  it.numberInSelf() :
                  it.numberInNeighbor()),
      dummy_(0.)
    {
    }
   
    //! copy constructor 
    ElementIntegrationPointList(const ElementIntegrationPointList& org)
      : quad_(org.quad_)
      , referenceGeometry_(org.referenceGeometry_)
      , elementGeometry_(org.elementGeometry_)
      , faceNumber_(org.faceNumber_)
      , dummy_(org.dummy_)
    {}
    
    //! The total number of integration points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith integration point.
    const CoordinateType& point(size_t i) const {
      dummy_ = referenceGeometry_.global(quad_.point(i));
      return dummy_;
    }

    //! Access to the ith integration point in local (codim-1 reference element)
    //! coordinates
    const LocalCoordinateType& localPoint(size_t i) const {
      return quad_.point(i);
    }

    //! A unique id per quadrature type.
    //! Quadratures are considered as distinct when they differ in the
    //! following points: geometry type, order, dimension and implementation.
    //! \note At the time of writing this, there is only one implementation
    //! per geometry type, order and dimension provided, but the concept is
    //! easily extendible beyond that.
    size_t id() const {
      return quad_.id();
    }

    //! The actual order of the quadrature.
    //! The actual order can be higher as the desired order when no 
    //! implementation for the desired order is found.
    int order() const {
      return quad_.order();
    }

    //! The geometry type the quadrature points belong to.
    GeometryType geometry() const {
      return quad_.geo();
    }

    //! The geometry type of the codim 0 reference element.
    GeometryType elementGeometry() const {
      return elementGeometry_;
    }

  protected:
    // return local face number 
    int faceNumber() const { return faceNumber_; }

    // return reference to integration point list 
    const IntegrationPointListType& quadImp() const { return quad_; }

  private:
   typedef typename IntersectionIterator::LocalGeometry ReferenceGeometry;

  protected:
    const IntegrationPointListType quad_;
    const ReferenceGeometry& referenceGeometry_;
    const GeometryType elementGeometry_;
    const int faceNumber_;

    mutable CoordinateType dummy_;
  };

} // end namespace Dune
#endif
