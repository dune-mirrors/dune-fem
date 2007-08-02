#ifndef DUNE_FEM_QUADRATURE_HH
#define DUNE_FEM_QUADRATURE_HH

//#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/misc/array.hh>

/*
#ifndef HAVE_ALBERTA_FOUND
#ifdef HAVE_ALBERTA
#if HAVE_ALBERTA 
#define HAVE_ALBERTA_FOUND
#warning "Use ALBERTA quadratures!"
#endif
#endif
#endif
*/

// use quadratures from dune-grid 
#define USE_DUNE_QUADRATURES

// do not use ALBERTA Quadratures at the moment 
#undef HAVE_ALBERTA_FOUND
#ifdef HAVE_ALBERTA_FOUND
// inlcude albertagrid.hh includes the needed alberta.h 
#include <dune/grid/albertagrid.hh>
#endif

// quadrature storage classes 
#include "quadprovider.hh"
#include "gausspoints.hh"
#include "prismpoints.hh"
#include "pyramidpoints.hh"

// include quadrature points 
#ifdef USE_DUNE_QUADRATURES
#include "dunequadratures.hh"
#else
#include "simplexpoints.hh"
#endif

namespace Dune
{

  // Forward declaration
  template <typename ct, int dim, template <class,int> class QuadratureTraits> 
  class QuadratureProvider;
  template <class ct, int dim>
  class QuadratureImp;
  
#ifndef USE_DUNE_QUADRATURES
  template <class ct, int dim>
  class SimplexQuadrature;
  template <class ct, int dim>
  class CubeQuadrature;
  template <class ct>
  class LineQuadrature;
  template <class ct>
  class TriangleQuadrature;
  template <class ct>
  class TetraQuadrature;
  template <class ct>
  class QuadrilateralQuadrature;
  template <class ct>
  class HexaQuadrature;
  template <class ct>
  class PrismQuadrature;
  template <class ct>
  class PyramidQuadrature;
#endif

  /*! \defgroup Quadrature Quadratures
   *
   *  In DUNE, quadratures are a set of quadrature points and corresponding
   *  weights.
   *
   *  \remark To get an overview of a quadrature, see Quadrature.
   *
   *  \remark The quadratures usually used are ElementQuadrature and
   *  CachingQuadrature.
   *
   *  \endgroup
   */



  /*! \class IntegrationPointListImp
   *  \ingroup Quadrature
   *  \brief Generic implementation of an IntegrationPointList
   *
   *  An integration point list is simply a list of points, given in local
   *  coordinates, i.e., coordinates within the reference element.
   *
   *  \note Integration point lists do not change over time. It can safely
   *        be assumed that they always return the same points in the same
   *        order.
   */
  template< typename FieldImp, int dim >
  class IntegrationPointListImp
  {
  public:
    //! field type
    typedef FieldImp FieldType;

  private:
    typedef IntegrationPointListImp< FieldType, dim > ThisType;

  public:
    //! Local coordinate type
    typedef FieldVector< FieldType, dim > CoordinateType;
    
    //! to be revised, look at caching quad 
    enum { codimension = 0 };
    
  private:
    // vector holding the coordinates for each point
    DynamicArray< CoordinateType > points_;
    //std :: vector< CoordinateType > points_;

    // identifier of the integration point list
    const size_t id_;
    
  protected:
    /*! \brief Constructor
     *
     *  The constructor simply creates an empty point list and stores the
     *  specified identifier.
     * 
     *  \note The constructors of derived classes should fill the integration
     *        point list via addIntegrationPoint.
     *
     *  \note The identifier of an integration point list must be globally
     *        unique. Even integration point lists for different dimensions
     *        must have different identifiers.
     *
     *  \param[in]  id  unique identifier of the integration point list
     *                  (provided by QuadratureProvider)
     */
    inline explicit IntegrationPointListImp( size_t id )
    : points_(),
      id_( id )
    {
    }
   
  private:
    // Copying is forbidden
    IntegrationPointListImp( const IntegrationPointListImp& );

  public:
    // Virtual destructor
    virtual ~IntegrationPointListImp ()
    {
    }

    /*! \brief obtain coordinates of i-th integration point
     *
     *  This method returns a reference to the coordinates of the i-th
     *  integration point for 0 <= i < nop(). The integration point is given
     *  in local coordinates, i.e., coordinates with respect to the reference
     *  element.
     * 
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns reference to i-th integration point
     */
    inline const CoordinateType &point ( size_t i ) const
    {
      assert( i < nop() );
      return points_[ i ];
    }

    /*! \brief obtain the number of integration points
     *
     *  \returns number of integration points within this list
     */
    size_t nop () const
    {
      return points_.size();
    }

    /*! \brief obtain the identifier of the integration point list
     * 
     *  The identifier of an integration point list must be globally unique.
     *  Even integration point lists for different dimensions must have
     *  different identifiers.
     *
     *  \note Quadratures are considered distinct if they differ in one of the
     *        following points: geometry type, order, dimension or implementation.
     * 
     *  \returns globally unique identifier of the integration point list
     */
    size_t id () const
    {
      return id_;
    }

    /*! \brief obtain order of the integration point list
     *
     *  The order of a quadrature is the maximal polynomial degree that is
     *  guaranteed to be integrated exactly by the quadrature.
     *
     *  In case of an integration point list, the definition of this value is
     *  left to the implementor.
     *
     *  \returns the order of the integration point list
     */
    virtual int order() const = 0;

    /*! \brief obtain GeometryType for this integration point list
     *
     *  Integration point lists are specified in local coordinates, i.e.,
     *  coordinates with respect to the reference element. Hence, each 
     *  integration point list is only valid for one type of geometry, i.e.,
     *  for one reference element. The type can be retrieved via this method.
     *
     *  \returns GeometryType for this integration point list
     */
    virtual GeometryType geometry () const = 0;

  protected:
    /*! \brief Adds an integration point to the list
     *
     *  This method allows derived classes to add integration points to the
     *  list. This mehtod should only be used within the constructor of the
     *  derived class.
     */
    void addIntegrationPoint( const CoordinateType &point )
    {
      //points_.push_back( point );
      points_.append( point );
    }
  };



  /*! \class QuadratureImp
   *  \ingroup Quadrature
   *  \brief Generic implementation of a Dune quadrature.
   *
   *  A Dune Quadrature is nothing but a list of integration points (see also
   *  IntegrationPointsListImp) and their respective weights.
   *
   *  \note Quadratures do not change over time. It can safely be assumed that
   *        they always return the same points in the same order.
   */
  template< typename FieldImp, int dim >
  class QuadratureImp
  : public IntegrationPointListImp< FieldImp, dim >
  {
  public:
    //! field type
    typedef FieldImp FieldType;
  
  private:
    typedef QuadratureImp< FieldType, dim > ThisType;
    typedef IntegrationPointListImp< FieldType, dim > BaseType;

  public:
    //! \copydoc Dune::IntegrationPointsListImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    // vector holding weights of each integration point 
    DynamicArray< FieldType > weights_;
    //std :: vector< FieldType > weights_;
 
  protected:
    /*! \brief Constructor
     *
     *  The constructor simply creates an empty quadrature and stores the
     *  specified identifier.
     * 
     *  \note The constructors of derived classes should fill the quadrature
     *        via addQuadraturePoint
     *
     *  \note The identifier of an integration point list must be globally
     *        unique. Even integration point lists for different dimensions
     *        must have different identifiers.
     *
     *  \param[in]  id  unique identifier of the quadrature (provided by
     *                  QuadratureProvider)
     */
    inline explicit QuadratureImp( size_t id )
    : BaseType( id ),
      weights_()
    {
    }
   
  private:
    // Copying is forbidden
    QuadratureImp ( const QuadratureImp& );

  public:
    virtual ~QuadratureImp ()
    {
    }
    
    /*! \brief obtain weight of i-th integration point
     *
     *  This method returns the weight of the i-th integration point for
     *  0 <= i < nop() within the quadrature.
     *
     *  \note The integration point can be obtained via the point() method.
     * 
     *  \note The quadrature weights sum up to the volume of the reference
     *        element.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns weight of the i-th integration point
     */
    const FieldType &weight ( size_t i ) const
    {
      return weights_[ i ];
    }

  private:
    // Disallow use of addIntegrationPoint for quadratures
    inline void addIntegrationPoint ( const CoordinateType &point )
    {
      BaseType :: addIntegrationPoint( point );
    }

  protected:
    /*! \brief Adds a point-weight pair to the quadrature
     *
     *  This method allows derived classes to add quadrature points (and their
     *  respective weights) to the list. This mehtod should only be used within
     *  the constructor of the derived class.
     */
    inline void addQuadraturePoint ( const CoordinateType &point,
                                     const FieldType weight )
    {
      addIntegrationPoint( point );
      //weights_.push_back( weight );
      weights_.append( weight );
    }
  };



#ifndef USE_DUNE_QUADRATURES
  /*! \class SimplexQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

#ifdef HAVE_ALBERTA_FOUND
    enum { maxOrder1 = 19, maxOrder2 = 17, maxOrder3 = 7 };
#else
    enum { maxOrder1 = 19, maxOrder2 = 12 , maxOrder3 = 5 };
#endif
          
  protected:
    int order_;
    
  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    SimplexQuadrature( const GeometryType& geometry, int order, size_t id );
    
    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, dim );
    }
   
    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures.
    static size_t maxOrder ()
    {
      if( dim == 1 )
        return maxOrder1;
      if( dim == 2 )
        return maxOrder2;
      if( dim == 3 )
        return maxOrder3;
      DUNE_THROW( NotImplemented, "SimplexQuadratures from dim > 3 not implemented." );
    }
  };
#endif



#ifndef USE_DUNE_QUADRATURES
  /*! \class CubeQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  protected:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    CubeQuadrature( const GeometryType &geometry, int order, size_t id );
    
    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, dim );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures.
    static size_t maxOrder ()
    { 
      return GaussPts :: highestOrder;
    }
  };
#endif
  


#ifndef USE_DUNE_QUADRATURES
  /*! \class LineQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  protected:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    LineQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 1 );
    }

    //! copydoc Dune::QuadratureImp::order
    virtual int order() const
    {
      return order_;
    }

    //! maximal order of available quadratures.
    static size_t maxOrder ()
    { 
      return GaussPts::highestOrder;
    }
  };
#endif
 


#ifndef USE_DUNE_QUADRATURES
  /*! \class TriangleQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    TriangleQuadrature ( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, 2 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures.
    static size_t maxOrder ()
    { 
#ifdef HAVE_ALBERTA_FOUND
      // highest order of Alberta quads 
      return 17; 
#else 
      // highest order of UG quads 
      return 12; 
#endif
    }
  };
#endif



#ifndef USE_DUNE_QUADRATURES
  /*! \class QuadrilateralQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    QuadrilateralQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 2 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures.
    static size_t maxOrder ()
    { 
      return GaussPts :: highestOrder;
    }
  };
#endif



#ifndef USE_DUNE_QUADRATURES
  /*! \class TetraQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;
    
  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    TetraQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: simplex, 3 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures
    static size_t maxOrder ()
    {
#ifdef HAVE_ALBERTA_FOUND
      // highest order of Alberta quads 
      return 7; 
#else 
      // highest order of UG quads 
      return 5; 
#endif
    }
  };
#endif



#ifndef USE_DUNE_QUADRATURES
  /*! \class HexaQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    HexaQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: cube, 3 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures
    static size_t maxOrder()
    { 
      return GaussPts::highestOrder;
    }
  };
#endif



#ifndef USE_DUNE_QUADRATURES
  /*! \class PrismQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType;
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    PrismQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: prism, 3 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures
    static size_t maxOrder () 
    {
      return PrismPoints :: highest_order;
    }
  };
#endif


  
#ifndef USE_DUNE_QUADRATURES
  /*! \class PyramidQuadrature
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
    //! \copydoc Dune::QuadratureImp::CoordinateType
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    int order_;

  public:
    /*! \brief constructor filling the list of points and weights
     *
     *  \param[in]  gemoetry  geometry type for which a quadrature is desired
     *  \param[in]  order     desired order (provided by the user)
     *  \param[in]  id        unique identifier (provided by QuadratureProvider)
     */
    PyramidQuadrature( const GeometryType &geometry, int order, size_t id );

    //! \copydoc Dune::QuadratureImp::geometry
    virtual GeometryType geometry () const
    {
      return GeometryType( GeometryType :: pyramid, 3 );
    }

    //! \copydoc Dune::QuadratureImp::order
    virtual int order () const
    {
      return order_;
    }

    //! maximal order of available quadratures
    static size_t maxOrder ()
    {
      return PyramidPoints :: highest_order;
    }
  };
#endif



  // \brief Allows injection of arbitrary points as quadrature points.
  // Useful to test some features of the quadrature framework in isolation
  // and with known input data. Each TestQuadrature object gets its own
  // unique id.
  template <class ct, int dim>
  class TestQuadrature : public QuadratureImp<ct, dim>
  {
  public:
    typedef FieldVector<ct, dim> CoordinateType;

    // dummy value
    enum { maxOrder_ = 10 };

  public:
    //! Constructor
    TestQuadrature(const GeometryType& geo, int order);

    //! Adds new quadrature point/weight pair
    void newQuadraturePoint(const CoordinateType& c, ct weight);

    //! Desired geometry
    virtual GeometryType geometry() const { return geo_; }

    //! Dummy order method
    virtual int order() const { return order_; }
    
    //! Dummy max order method
    static size_t maxOrder() { return maxOrder_; }
    
  private:
    GeometryType geo_;
    int order_;
  };

  

  // default defines for used quadratures 
  template< typename FieldType, int dim >
  struct DefaultQuadratureTraits
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, dim > CubeQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, dim > CubeQuadratureType; 
#endif
    typedef QuadratureImp< FieldType, dim > IntegrationPointListType;
  }; 



  // quadratures for points 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 0 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 0 > PointQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 0 > PointQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 0 > IntegrationPointListType;
  };
 


  // quadratures for lines 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 1 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 1 > LineQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 1 > LineQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 1 > IntegrationPointListType;
  };
 


  // quadratures for simplex and cubes 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 2 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 2 > SimplexQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 2 > CubeQuadratureType;
#else
    typedef CubeQuadrature< FieldType, 2 > CubeQuadratureType; 
    typedef SimplexQuadrature< FieldType, 2 > SimplexQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 2 > IntegrationPointListType;
  };


  
  // quadratures for simplex, cubes, prisms, and pyramids
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType , 3 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 3 > SimplexQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 3 > CubeQuadratureType;

    typedef QuadratureRulesFactory< FieldType, 3 > PrismQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 3 > PyramidQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 3 > CubeQuadratureType;     
    typedef SimplexQuadrature< FieldType, 3 > SimplexQuadratureType;     

    typedef PrismQuadrature< FieldType > PrismQuadratureType;
    typedef PyramidQuadrature< FieldType > PyramidQuadratureType;
#endif

    typedef QuadratureImp< FieldType, 3 > IntegrationPointListType;
  };



  /*! \class IntegrationPointList
   *  \ingroup Quadrature
   *  \brief actual interface class for integration point lists
   *
   *  IntegrationPointList is a proxy for the actual implementations of the
   *  integration point lists. During construction, the IntegrationPointList
   *  object is configured with an appropriate implementation object from the
   *  QuadratureProvider (monostate pattern).
   *
   *  The design goal is minimization of construction time. The actual
   *  implementation can be created once and reused whenever it is needed.
   *  Moreover, this layout insulates the user from all initialization and
   *  storage stuff.
   *
   *  \note The difference between integration point lists and quadratures is
   *        that quadratures have weights.
   */
  template< typename FieldImp, int dim,
            template< class, int > class IntegrationTraits >
  class IntegrationPointList 
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

  private:
    typedef IntegrationPointList< FieldType, dimension, IntegrationTraits > ThisType;

    typedef IntegrationTraits< FieldType, dimension > Traits;

    typedef QuadratureProvider< FieldType, dimension, IntegrationTraits >
      QuadratureProviderType;

  public:
    //! type of integration point list implementation 
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;
    
    //! type of coordinate
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  protected:
    const IntegrationPointListType &ipList_;

  public:
    /*! \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  order     order of the requested quadrature
     */
    inline IntegrationPointList( const GeometryType &geometry, int order )
    : ipList_( QuadratureProviderType :: getQuadrature( geometry, order ) )
    {
    }

    /*! \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for testing purposes.
     *
     *  \param[in]  ipList  implementation of the integration point list
     */
    inline IntegrationPointList(const IntegrationPointListType &ipList )
    : ipList_( ipList )
    {
    }

    /*! \brief copy constructor
     *  
     *  \param[in]  org  integration point list to be copied
     */ 
    inline IntegrationPointList( const IntegrationPointList &org )
    : ipList_( org.ipList_ )
    {
    }

    /*! \brief obtain a reference the actual implementation
     * 
     *  \returns a reference to the implementation of this integration point
     *           list
     */
    const IntegrationPointListType &ipList () const
    {
      return ipList_;
    }

    /*! \brief obtain the number of integration points
     *
     *  \returns number of integration points within this list
     */
    int nop () const
    {
      return ipList_.nop();
    }
    
    /*! \brief obtain coordinates of i-th integration point
     *
     *  This method returns a reference to the coordinates of the i-th
     *  integration point for 0 <= i < nop(). The integration point is given
     *  in local coordinates, i.e., coordinates with respect to the reference
     *  element.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns reference to i-th integration point
     */
    const CoordinateType &point ( size_t i ) const
    {
      return ipList_.point( i );
    }

    /*! \brief obtain the identifier of the integration point list
     * 
     *  The identifier of an integration point list must be globally unique.
     *  Even integration point lists for different dimensions must have
     *  different identifiers.
     *
     *  \note Quadratures are considered distinct if they differ in one of the
     *        following points: geometry type, order, dimension or implementation.
     *
     *  \returns globally unique identifier of the integration point list
     */
    size_t id () const
    {
      return ipList_.id();
    }
    
    /*! \brief obtain order of the integration point list
     *
     *  The order of a quadrature is the maximal polynomial degree that is
     *  guaranteed to be integrated exactly by the quadrature.
     *
     *  In case of an integration point list, the definition of this value is
     *  left to the implementor.
     *  
     *  \note Calling this method yields a virtual function call, so do not
     *        call this method unnecessarily.
     *
     *  \returns the order of the integration point list
     */
    int order () const
    {
      return ipList_.order();
    }
    
    /*! \brief obtain GeometryType for this integration point list
     *
     *  Integration point lists are specified in local coordinates, i.e.,
     *  coordinates with respect to the reference element. Hence, each 
     *  integration point list is only valid for one type of geometry, i.e.,
     *  for one reference element. The type can be retrieved via this method.
     *
     *  \note Calling this method yields a virtual function call, so do not
     *        call this method unnecessarily.
     *
     *  \returns GeometryType for this integration point list
     */
    GeometryType geometry () const
    {
      return ipList_.geometry();
    }
  };



  /*! \class Quadrature
   *  \ingroup Quadrature
   *  \brief actual interface class for quadratures
   *
   *  IntegrationPointList is a proxy for the actual implementations of the
   *  integration point lists. During construction, the IntegrationPointList
   *  object is configured with an appropriate implementation object from the
   *  QuadratureProvider (monostate pattern).
   *
   *  The design goal is minimization of construction time. The actual
   *  implementation can be created once and reused whenever it is needed.
   *  Moreover, this layout insulates the user from all initialization and
   *  storage stuff.
   *
   *  \note The difference between integration point lists and quadratures is
   *        that quadratures have weights.
   */
  template< class FieldImp, int dim,
            template< class, int > class QuadratureTraits = DefaultQuadratureTraits >
  class Quadrature
  : public IntegrationPointList< FieldImp, dim, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

  private:
    typedef Quadrature< FieldType, dimension, QuadratureTraits > ThisType;
    typedef IntegrationPointList< FieldType, dimension, QuadratureTraits > BaseType;

    typedef QuadratureTraits< FieldType, dimension > Traits;
    
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits >
      QuadratureProviderType;

  public:
    using BaseType :: ipList;

    //! type of the implementation (this must actually be a quadrature implementation)
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;

    //! type of local coordinate vectors
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  public:
    /*! \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  order     order of the requested quadrature
     */
    inline Quadrature( const GeometryType &geometry, int order )
    : BaseType( geometry, order )
    {
    }

    /*! \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for testing purposes.
     *
     *  \param[in]  ipList  implementation of the integration point list
     */
    inline explicit Quadrature( const IntegrationPointListType& ipList )
    : BaseType( ipList )
    {
    }

    /*! \brief copy constructor
     *  
     *  \param[in]  org  quadrature to be copied
     */ 
   //! Copy constructor
    inline Quadrature( const Quadrature &org )
    : BaseType( org )
    {
    }

    /*! \brief obtain weight of i-th integration point
     *
     *  This method returns the weight of the i-th integration point for
     *  0 <= i < nop() within the quadrature.
     *
     *  \note The integration point can be obtained via the point() method.
     *
     *  \note The quadrature weights sum up to the volume of the reference
     *        element.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns weight of the i-th integration point
     */
    const FieldType &weight( size_t i ) const
    {
      return ipList().weight( i );
    }
  };

} // end namespace Dune

#include "quadrature.cc"

#undef USE_DUNE_QUADRATURES
#undef HAVE_ALBERTA_FOUND

#endif
