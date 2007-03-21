#ifndef DUNE_LAGRANGEPOINTS_HH
#define DUNE_LAGRANGEPOINTS_HH

#include <dune/grid/common/geometry.hh>

#include <dune/fem/quadrature/quadrature.hh>

namespace Dune {
   
  /** \brief Actual implementation of LagrangeQuadrature for Geometry type Line
   */
  template< class ct >
  class LagrangeLineQuadrature : public QuadratureImp< ct, 1 >
  {
  public:
    //! Remember the field type
    typedef ct ctype;
    
    //! Type for coordinates
    typedef FieldVector< ctype, 1 > CoordinateType;

  private:
    typedef LagrangeLineQuadrature< ctype > ThisType;

    const int order_;
    
  public:
    LagrangeLineQuadrature( int order, size_t id );

    virtual GeometryType geometry() const
    {
      return GeometryType( GeometryType :: cube, 1 );
    }
    
    /** \brief Order of the quadrature
     * 
     *  Returns the order of the associated Lagrange polynom (as requested in
     *  the constructor).
     */
    virtual int order() const
    {
      return order_;
    }

    /** \brief Maximal order for this quadrature
     * 
     * Returns the maximal order, this quadrature can be constructed with. For
     * now, this is 2 (i.e. quadratic polynoms).
     */
    static size_t maxOrder()
    {
      return 2;
    }
  };



  /** \brief Actual implementation of LagrangeQuadrature for Geometry type
   *  Triangle
   */
  template< class ct >
  class LagrangeTriangleQuadrature : public QuadratureImp< ct, 2 >
  {
  public:
    //! Remember the field type
    typedef ct ctype;
    
    //! Type for coordinates
    typedef FieldVector< ctype, 2 > CoordinateType;

  private:
    typedef LagrangeTriangleQuadrature< ctype > ThisType;

    const int order_;
    
  public:
    LagrangeTriangleQuadrature( int order, size_t id );

    virtual GeometryType geometry() const
    {
      return GeometryType( GeometryType :: simplex, 2 );
    }
    
    /** \brief Order of the quadrature
     * 
     *  Returns the order of the associated Lagrange polynom (as requested in
     *  the constructor).
     */
    virtual int order() const
    {
      return order_;
    }

    /** \brief Maximal order for this quadrature
     * 
     * Returns the maximal order, this quadrature can be constructed with. For
     * now, this is 2 (i.e. quadratic polynoms).
     */
    static size_t maxOrder()
    {
      return 2;
    }
  };

  

  /** \brief Actual implementation of LagrangeQuadrature for Geometry type
   *  Quadrilateral
   */
  template< class ct >
  class LagrangeQuadrilateralQuadrature : public QuadratureImp< ct, 2 >
  {
  public:
    //! Remember the field type
    typedef ct ctype;
    
    //! Type for coordinates
    typedef FieldVector< ctype, 2 > CoordinateType;

  private:
    typedef LagrangeQuadrilateralQuadrature< ctype > ThisType;

    const int order_;
    
  public:
    LagrangeQuadrilateralQuadrature( int order, size_t id );

    virtual GeometryType geometry() const
    {
      return GeometryType( GeometryType :: cube, 2 );
    }
    
    /** \brief Order of the quadrature
     * 
     *  Returns the order of the associated Lagrange polynom (as requested in
     *  the constructor).
     */
    virtual int order() const
    {
      return order_;
    }
    
    /** \brief Maximal order for this quadrature
     * 
     * Returns the maximal order, this quadrature can be constructed with. For
     * now, this is 2 (i.e. quadratic polynoms).
     */
    static size_t maxOrder()
    {
      return 2;
    }
  };



  //! Traits for Lagrange quadratures
  template< class ct, int dim >
  struct LagrangeQuadratureTraits
  {
    typedef CompileTimeChecker< false > Only_implementations_for_dim_1_2_exist;
  };

  /*
  template< class ct >
  struct LagrangeQuadratureTraits< ct, 0 >
  {
    typedef CubeQuadrature< ct, 0 > PointQuadratureType;
  };
  */

  //! Lagrange quadratures for lines
  template< class ct >
  struct LagrangeQuadratureTraits< ct, 1 >
  {
    typedef LagrangeLineQuadrature< ct > LineQuadratureType;
  };

  //! Lagrange quadratures for simplex and cubes
  template< class ct >
  struct LagrangeQuadratureTraits< ct, 2 >
  {
    typedef LagrangeQuadrilateralQuadrature< ct > CubeQuadratureType;
    typedef LagrangeTriangleQuadrature< ct >      SimplexQuadratureType;
  };

  /*
  template< class ct >
  struct LagrangeQuadratureTraits< ct, 3 >
  {
    typedef CubeQuadrature< ct, 3 >    CubeQuadratureType;
    typedef SimplexQuadrature< ct, 3 > SimplexQuadratureType;

    typedef PrismQuadrature< ct >      PrismQuadratureType;
    typedef PyramidQuadrature< ct >    PyramidQuadratureType;
  };
  */



  /** \brief Lagrange Quadrature
   *
   * Sometimes it is necessary to know the coordinate of a Lagrange point.
   * This quadrature enumerates these points (in the same order as the local
   * DoFs of a Lagrange polynomial).
   *
   * \note The order of this quadrature corresponds to the order of the
   * Lagrange polynoms, which is less than the order of the quadrature
   * (when used as a quadrature for integration).
   *
   * \note Only LagrangeQuadrature< GridPartImp, 0 > is implemented.
   */
  template< typename GridPartImp, int codim >
  class LagrangeQuadrature
  {
    typedef CompileTimeChecker< false > Only_implementation_for_codim_0_exists;
  };

  /** \brief Lagrange Quadrature (for codimension 0)
   *
   * Sometimes it is necessary to know the coordinate of a Lagrange point.
   * This quadrature enumerates these points (in the same order as the local
   * DoFs of a Lagrange polynomial).
   *
   * \note The order of this quadrature corresponds to the order of the
   * Lagrange polynoms, which is less than the order of the quadrature
   * (when used as a quadrature for integration).
   */
  template< typename GridPartImp >
  class LagrangeQuadrature< GridPartImp, 0 >
  {
  public:
    //! Remember the GridPart
    typedef typename GridPartImp :: GridType GridPartType;
    
     //! The type for reals (usually double)
    typedef typename GridPartType :: ctype ctype;
    
    //! Dimension of the world
    enum { dimension = GridPartType :: dimension };
    
    //! Codimension is zero by definition
    enum { codimension = 0 };

    enum Side { INSIDE, OUTSIDE };
    
    //! Type of the codim-0 entity
    typedef typename GridPartType :: template Codim< 0 > :: Entity EntityType;

  private:
    typedef Quadrature< ctype, dimension, LagrangeQuadratureTraits >
      QuadratureType;

    typedef LagrangeQuadrature< GridType, 0 > ThisType;
  
  public:
    //! Type for coordinates in the codim-0 reference element
    typedef typename QuadratureType :: CoordinateType CoordinateType;

  private:
    QuadratureType quad_;
    
  public:
    //! Constructor
    //! \param entity Entity the quadrature lives on (respectively on its reference element).
    //! \param order Desired minimal order of the quadrature.
    LagrangeQuadrature( const EntityType &entity, int order )
      : quad_( entity.geometry().type(), order )
    {
    }

    //! Number of Lagrange points
    int nop() const
    {
      return quad_.nop();
    }

    //! Access the i-th Lagrange point
    const CoordinateType& point( size_t i ) const
    {
      return quad_.point( i );
    }

    //! Access the i-th Lagrange point
    const CoordinateType& localPoint( size_t i ) const
    {
      return quad_.point( i );
    }

    /** \brief The weight of the i-th quadrature point
     *
     *  The quadrature weights sum up to the volume of the respective reference
     *  element.
     *
     *  \note All Lagrange points have the same weight.
     */
    const ctype& weight( size_t i ) const
    {
      return quad_.weight( i ) ;
    }

    //! A unique id per quadrature type.
    //! Quadratures are considered as distinct when they differ in the
    //! following points: geometry type, order, dimension and implementation.
    //! \note At the time of writing this, there is only one implementation
    //! per geometry type, order and dimension provided, but the concept is
    //! easily extendible beyond that.
    size_t id() const
    {
      return quad_.id();
    }

    /** \brief Order of the quadrature
     * 
     *  Returns the order of the associated Lagrange polynom (as requested in
     *  the constructor).
     */
    int order() const
    {
      return quad_.order();
    }

    //! Type of geometry the quadrature points belong to
    GeometryType geometry() const
    {
      return quad_.geometry();
    }

    //! Geometry type of the associated element of codimension 0
    GeometryType elementGeometry() const
    {
      return quad_.geometry();
    }

    //! Returns quadraturePoint to behave like a caching qaudrature (but
    //! without caching)
    //! \note This only works for codimension 0
    size_t cachingPoint( size_t quadraturePoint ) const
    {
      return quadraturePoint;
    }

  protected:
    const QuadratureType& quadImp() const
    {
      return quad_;
    }
  };

}

#include "lagrangepoints.cc"
#endif // DUNE_LEGRANGEPOINTS_HH
