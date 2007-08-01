#ifndef DUNE_ELEMENTQUADRATURE_HH
#define DUNE_ELEMENTQUADRATURE_HH

#include "quadrature.hh"
#include "elementpointlist.hh"

namespace Dune
{
  
  /*! \class ElementQuadrature
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
   *  corresponding to the geometry to the codim-0 reference element. Moreover,
   *  the transformation is done only once (when the quadrature is created).
   *  
   *  To achieve this goal, an element quadrature depends stronger on the
   *  context in which it is used. For example, for each face within a
   *  tetrahedron (though they are all the same) we need a different
   *  ElementQuadrature, since the coordinates of the quadrature points
   *  with respect to the codim-0 entity differ for each face.
   *
   *  \note Actually, codim-1 element quadratures depend on the intersection.
   */
  template< typename GridPartImp, int codim >
  class ElementQuadrature
  {
    typedef CompileTimeChecker< false >
      __Only_implementations_for_codim_0_and_1_exist__;
  };



  template< class GridPartImp, int codim >
  struct ElementQuadratureTraits
  {
    // type of single coordinate
    typedef typename GridPartImp :: GridType :: ctype ctype;

    // dimension of quadrature 
    enum { dimension = GridPartImp :: GridType :: dimension };

    // codimension of quadrature
    enum { codimension = codim };

    // type of used integration point list
    typedef Quadrature< ctype, dimension-codim, DefaultQuadratureTraits > IntegrationPointListType;
    
    // type of local coordinate (with respect to the codim-0 entity)
    typedef typename Quadrature< ctype, dimension, DefaultQuadratureTraits > :: CoordinateType
      CoordinateType; 
  };

  

  // \brief Element quadrature on codim-0 entities.
  // For codim-0 element quadratures, there is no additional information
  // from the context needes, in consequence, the quadrature behaves like
  // a generic quadrature class, independent from the situation in the grid.
  template< typename GridPartImp >
  class ElementQuadrature< GridPartImp, 0 >
  : public ElementIntegrationPointList
    < GridPartImp, 0, ElementQuadratureTraits< GridPartImp, 0 > >
  {
  public:
    //! type of the grid partition
    typedef GridPartImp GridPartType;

    //! codimension of the element quadrature
    enum { codimension = 0 };

  private:
    typedef ElementQuadratureTraits< GridPartType, codimension > IntegrationTraits;
    
    typedef ElementQuadrature< GridPartType, codimension > ThisType;
    typedef ElementIntegrationPointList< GridPartType, codimension, IntegrationTraits >
      BaseType;

    using BaseType :: quadImp;
  
  public:
    //! type of the grid
    typedef typename GridPartType :: GridType GridType;

    //! dimension of the world
    enum { dimension = GridType :: dimension };
    
    //! type for reals (usually double)
    typedef typename GridType :: ctype RealType;
    
  public:
    //! type for coordinates in the codim-0 reference element 
    typedef typename IntegrationTraits :: CoordinateType CoordinateType;
    
    //! type of the codim-0 entity
    typedef typename GridType :: template Codim< 0 > :: Entity Entity;

  public:
    /*! \brief copy constructor
     *
     *  \param[in]  org  element quadrature to copy
     */
    ElementQuadrature( const ElementQuadrature& org )
    : BaseType( org )
    {
    }
   
    /*! \brief constructor
     *
     *  \param[in]  entity  entity, on whose reference element the quadratre
     *                      lives
     *  \param[in]  order   desired minimal order of the quadrature
     */
    ElementQuadrature( const Entity &entity, int order )
    : BaseType( entity.geometry().type(), order )
    {
    }
   
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



  // \brief Element quadrature on codim-1 entities.
  // For codimension 1, the quadrature needs information about the 
  // intersection. Plus, the user must decide if the quadrature shall live
  // on the reference element of the outside or inside element of the 
  // intersection.
  template< class GridPartImp >
  class ElementQuadrature< GridPartImp, 1 >
  : public ElementIntegrationPointList
    < GridPartImp, 1, ElementQuadratureTraits< GridPartImp, 1 > >
  {
  public:
    //! type of the grid partition
    typedef GridPartImp GridPartType;

    //! codimension of the quadrature
    enum { codimension = 1 };

  private:
    typedef ElementQuadratureTraits< GridPartType, codimension > IntegrationTraits;
    
    typedef ElementQuadrature< GridPartType, codimension > ThisType;
    typedef ElementIntegrationPointList< GridPartType, 1, IntegrationTraits >
      BaseType;

    using BaseType :: quadImp;
    
  public:
    //! type of the grid
    typedef typename GridPartType :: GridType GridType;
  
    //! dimension of the world
    enum { dimension = GridType :: dimension };
    
    //! type for reals (usually double)
    typedef typename GridType :: ctype RealType;

    //! type of the intersection iterator
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

    // just for compatibility (shoud be removed)
    typedef IntersectionIteratorType IntersectionIterator;

    //! type of coordinates in codim-0 reference element
    typedef typename IntegrationTraits :: CoordinateType CoordinateType;

    //! type of coordinate in codim-1 reference element
    typedef typename IntegrationTraits :: IntegrationPointListType :: CoordinateType
      LocalCoordinateType;

    //! type of quadrature for use on non-conforming intersections 
    typedef ThisType NonConformingQuadratureType;
    
  public:
    /*! \brief constructor
     *
     *  \param[in]  gridPart      grid partition (a dummy here)
     *  \param[in]  intersection  intersection iterator
     *  \param[in]  order         desired order of the quadrature
     *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
     *                            which the ElementQuadrature shall be created
     */
    ElementQuadrature ( const GridPartType &gridPart,
                        const IntersectionIteratorType &intersection, 
                        int order, 
                        typename BaseType :: Side side )
    : BaseType( gridPart, intersection, order, side )
    {
    }

    /*! \brief constructor
     *
     *  \param[in]  intersection  intersection iterator
     *  \param[in]  order         desired order of the quadrature
     *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
     *                            which the ElementQuadrature shall be created
     */
    ElementQuadrature( const IntersectionIterator &intersection,
                       int order,
                       typename BaseType :: Side side )
    : BaseType( intersection, order, side )
    {
    }
    
    /*! \brief copy constructor
     *
     *  \param[in]  org  element quadrature to copy
     */
    ElementQuadrature( const ElementQuadrature &org )
    : BaseType( org )
    {
    }
    
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

} // end namespace Dune
#endif
