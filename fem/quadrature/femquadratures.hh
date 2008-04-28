#ifndef DUNE_FEM_FEMQUADRATURES_HH
#define DUNE_FEM_FEMQUADRATURES_HH

#include <dune/common/geometrytype.hh>

#include <dune/fem/quadrature/quadratureimp.hh>

// quadrature storage classes 
#include "gausspoints.hh"
#include "pyramidpoints.hh"
#include "simplexpoints.hh"

namespace Dune
{

  struct SimplexMaxOrder 
  {
    // uses implementation from parDG
    enum { maxOrder1 = 39, maxOrder2 = 13, maxOrder3 = 12 };
  };
          
  /*  \class SimplexQuadrature
   *  \ingroup Quadrature
   *  \brief generic quadrature class for simplices
   *  
   *  SimplexQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The UG quadrature rules are used here. 
   */
  template< class FieldImp, int dim >
  class SimplexQuadrature
  : public QuadratureImp< FieldImp, dim >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef SimplexQuadrature< FieldType, dim > ThisType;
    typedef QuadratureImp< FieldType, dim > BaseType;

  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  protected:
    int order_;
    
  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    SimplexQuadrature( const GeometryType& geometry, int order, size_t id );
    
    /** \copydoc Dune::QuadratureImp::geometry
     */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, dim );
    }
   
    /** \copydoc Dune::QuadratureImp::order
     */
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures
    static size_t maxOrder ()
    {
      if( dim == 1 )
        return SimplexMaxOrder::maxOrder1;
      if( dim == 2 )
        return SimplexMaxOrder::maxOrder2;
      if( dim == 3 )
        return SimplexMaxOrder::maxOrder3;
      DUNE_THROW( NotImplemented, "SimplexQuadratures from dim > 3 not implemented." );
    }
  };



  /*  \class CubeQuadrature
   *  \ingroup Quadrature
   *  \brief generic quadrature class for cubes
   *  
   *  CubeQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The quadrature uses the 1d gauss points (and their tensorial
   *        product) as quadrature points
   */
  template< class FieldImp, int dim >
  class CubeQuadrature
  : public QuadratureImp< FieldImp, dim >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef CubeQuadrature< FieldType, dim > ThisType;
    typedef QuadratureImp< FieldType, dim > BaseType;
  
  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  protected:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    CubeQuadrature( const GeometryType &geometry, int order, size_t id );
    
    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, dim );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder ()
    { 
      return GaussPts :: highestOrder;
    }
  };
  


  /*  \class LineQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for lines
   *  
   *  LineQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note This class is redundant as CubeQuadrature can be used instead
   */
  template< class FieldImp >
  class LineQuadrature
  : public QuadratureImp< FieldImp, 1 > 
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef LineQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 1 > BaseType;
    
  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  protected:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    LineQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 1 );
    }

    /** copydoc Dune::QuadratureImp::order */
    virtual int order() const
    {
      return order_;
    }

    /** \brief  maximal order of available quadratures */
    static size_t maxOrder ()
    { 
      return GaussPts::highestOrder;
    }
  };
 


  /*  \class TriangleQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for triangles
   *  
   *  TriangleQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The UG quadrature rules are used here. 
   *  
   *  \note This class is redundant as SimplexQuadrature can be used instead.
   */
  template< class FieldImp >
  class TriangleQuadrature
  : public QuadratureImp< FieldImp, 2 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef TriangleQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 2 > BaseType;
    
  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    TriangleQuadrature ( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, 2 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder ()
    { 
      return SimplexMaxOrder::maxOrder2;
    }
  };



  /*  \class QuadrilateralQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for quadrilaterals
   *  
   *  QuadrilateralQuadrature implements the geometry-specific part of the
   *  quadrature and initialises the vector quadrature points and weights.
   *  
   *  \note The quadrature uses tensorial products of the 1d gauss points
   *        as quadrature points.
   *
   *  \note This class is redundant as CubeQuadrature can be used instead.
   */
  template< class FieldImp >
  class QuadrilateralQuadrature
  : public QuadratureImp< FieldImp, 2 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef QuadrilateralQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 2 > BaseType;
    
  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    QuadrilateralQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 2 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder ()
    { 
      return GaussPts :: highestOrder;
    }
  };



  /*  \class TetraQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for tetrahedra
   *  
   *  TetraQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The UG quadrature rules are used here. 
   *  
   *  \note This class is redundant as SimplexQuadrature can be used instead.
   */
  template< class FieldImp >
  class TetraQuadrature
  : public QuadratureImp< FieldImp, 3 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef TetraQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 3 > BaseType;

  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    TetraQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, 3 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder ()
    {
      return SimplexMaxOrder::maxOrder3;
    }
  };



  /*  \class HexaQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for hexahedra
   *  
   *  HexaQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The quadrature uses tensorial products of the 1d gauss points
   *        as quadrature points.
   *
   *  \note This class is redundant as CubeQuadrature can be used instead.
   */
  template< class FieldImp >
  class HexaQuadrature
  : public QuadratureImp< FieldImp, 3 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef HexaQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 3 > BaseType;

  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    HexaQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 3 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder()
    { 
      return GaussPts::highestOrder;
    }
  };



  /*  \class PrismQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for prisms
   *  
   *  PrismQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The HD stuff is used here, but needs some rework since only one
   *        rule is provided.
   */
  template< class FieldImp >
  class PrismQuadrature
  : public QuadratureImp< FieldImp, 3 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef PrismQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 3 > BaseType;
    
  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    PrismQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: prism, 3 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder () 
    {
      return SimplexMaxOrder::maxOrder2;
    }
  };


  
  /*  \class PyramidQuadrature
   *  \ingroup Quadrature
   *  \brief quadrature class for pyramids
   *  
   *  PyramidQuadrature implements the geometry-specific part of the quadrature
   *  and initialises the vector quadrature points and weights.
   *  
   *  \note The HD stuff is used here, but needs some rework since only one
   *        rule is provided.
   */
  template< class FieldImp >
  class PyramidQuadrature
  : public QuadratureImp< FieldImp, 3 >
  {
  public:
    typedef FieldImp FieldType;

  private:
    typedef PyramidQuadrature< FieldType > ThisType;
    typedef QuadratureImp< FieldType, 3 > BaseType;

  public:
    /** \copydoc Dune::QuadratureImp::CoordinateType */
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /** \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    PyramidQuadrature( const GeometryType &geometry, int order, size_t id );

    /** \copydoc Dune::QuadratureImp::geometry */
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: pyramid, 3 );
    }

    /** \copydoc Dune::QuadratureImp::order */
    virtual int order () const
    {
      return order_;
    }

    /** \brief maximal order of available quadratures */
    static size_t maxOrder ()
    {
      return PyramidPoints :: highest_order;
    }
  };

} // end namespace Dune

#include "femquadratures_inline.hh"

#endif
