#ifndef DUNE_FEMQUADRATURE_HH
#define DUNE_FEMQUADRATURE_HH

#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

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

  /*! \class IntegrationPointListImp
   *  \brief Generic implementation of an IntegrationPointList
   *
   *  An integration point list is simply a list of points, given in local
   *  coordinates, i.e., coordinates within the reference element.
   *
   *  \note Integration point lists do not change over time. It can safely
   *        be assume that they always return the same points in the same
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
    std :: vector< CoordinateType > points_;

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
     */
    inline IntegrationPointListImp( size_t id )
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
     *  \note The identifier of an integration point list must be globally
     *        unique. Even integration point lists for different dimensions
     *        must have different identifiers.
     * 
     *  \returns globally unique identifier of the integration point list
     */
    size_t id () const
    {
      return id_;
    }

    // This method belongs into quadrature!
    //! Maximal degree of polynomial that gets integrated exactly by the
    //! quadrature.
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
      points_.push_back( point );
    }
  };



  /*! \class QuadratureImp
   *  \brief Generic implementation of a Dune quadrature.
   *
   *  A Dune Quadrature is nothing but a list of integration points (see also
   *  IntegrationPointsListImp) and their respective weights.
   *
   *  \note Quadratures do not change over time. It can safely be assume that
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
    //! Local coordinate type
    typedef typename BaseType :: CoordinateType CoordinateType;

  private:
    // vector holding weights of each integration point 
    std :: vector< FieldType > weights_;
 
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
     *  \param[in]  id  unique identifier of the quadrature
     */
    inline QuadratureImp( size_t id )
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
      weights_.push_back( weight );
    }
  };



#ifndef USE_DUNE_QUADRATURES
  //! A generic quadrature class for simplices.
  //! The UG quadrature rules are used here. SimplexQuadrature implements
  //! the geometry-specific part of the quadrature and initialises the vector
  //! of quadrature points and weights.
  template <class ct, int dim>
  class SimplexQuadrature : public QuadratureImp<ct, dim>
  {
  public:
    typedef FieldVector<ct, dim> CoordinateType;

    enum { 
#ifdef HAVE_ALBERTA_FOUND
      maxOrder2 = 17 , maxOrder3 = 7 
#else 
      maxOrder2 = 12 , maxOrder3 = 5 
#endif
    };
          
    enum { maxOrder_  = (dim == 1) ? 19 : ((dim == 2) ? maxOrder2 : maxOrder3 ) };
  public:
    //! Constructor
    //! \param order The desired order (provided by the user)
    //! \param id A unique id (provided by QuadratureProvider)
    inline
    SimplexQuadrature(const GeometryType&, int order, size_t id);
    
    //! The geometry type is... simplex!
    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::simplex,dim);
    }
    
    //! Returns the effective order of the quadrature
    //! (can be higher than the desired)
    virtual int order() const {
      return order_;
    }

    //! The maximal order of simplex quadratures. This is to be understood
    //! as an upper bound...
    static size_t maxOrder() { return maxOrder_; }

  private:
    int order_;
  };

  //! A generic quadrature for cubes
  //! This quadrature uses the 1d gauss points (and the tensorial product
  //! thereof) as quadrature points.
  template <class ct, int dim>
  class CubeQuadrature : public QuadratureImp<ct, dim>
  {
  public:
    typedef FieldVector<ct, dim> CoordinateType;

  public:
    //! Constructor
    //! \param order The desired order (provided by the user)
    //! \param id A unique id (provided by QuadratureProvider)
    CubeQuadrature(const GeometryType&, int order, size_t id);
    
    //! The geometry type is... cube!
    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::cube,dim);
    }

    //! Returns the effective order of the quadrature
    //! (can be higher than the desired)
    virtual int order() const {
      return order_;
    }

    //! The maximal order of simplex quadratures.
    static size_t maxOrder() { return GaussPts::highestOrder; }

  private:
    int order_;
  };

  //! A quadrature class for lines
  //! \note This class is redundant as CubeQuadrature can be used instead
  template <class ct>
  class LineQuadrature : public QuadratureImp<ct, 1> 
  {
  public:
    typedef FieldVector<ct, 1> CoordinateType;

  public:
    LineQuadrature(const GeometryType&, int order, size_t id);

    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::cube,1);
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPts::highestOrder; }

  private:
    int order_;
  };

  //! A quadrature class for triangles
  //! \note This class is redundant as SimplexQuadrature can be used instead.
  template <class ct>
  class TriangleQuadrature : public QuadratureImp<ct, 2>
  {
  public:
    typedef FieldVector<ct, 2> CoordinateType;

  public:
    TriangleQuadrature(const GeometryType&, int order, size_t id);

    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::simplex,2);
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { 
#ifdef HAVE_ALBERTA_FOUND
      // highest order of Alberta quads 
      return 17; 
#else 
      // highest order of UG quads 
      return 12; 
#endif
    }

  private:
    int order_;
  };

  //! A quadrature class for quadrilaterals
  //! \note This class is redundant as CubeQuadrature can be used instead.
  template <class ct>
  class QuadrilateralQuadrature : public QuadratureImp<ct, 2>
  {
  public:
    typedef FieldVector<ct, 2> CoordinateType;

  public:
    QuadrilateralQuadrature(const GeometryType&, int order, size_t id);

    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::cube,2);
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPts::highestOrder; }

  private:
    int order_;
  };

  //! A quadrature class for tetrahedra
  //! \note This class is redundant as SimplexQuadrature can be used instead.
  template <class ct>
  class TetraQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    TetraQuadrature(const GeometryType&, int order, size_t id);

    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::simplex,3);
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { 
#ifdef HAVE_ALBERTA_FOUND
      // highest order of Alberta quads 
      return 7; 
#else 
      // highest order of UG quads 
      return 5; 
#endif
    }

  private:
    int order_;
  };

  //! A quadrature class for hexahedra
  //! \note This class is redundant as CubeQuadrature can be used instead.
  template <class ct>
  class HexaQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    HexaQuadrature(const GeometryType&, int order, size_t id);

    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::cube,3);
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPts::highestOrder; }

  private:
    int order_;
  };

  //! A quadrature class for prisms
  //! The HD stuff is used here, but needs some rework since only one rule
  //! is provided. But since nobody here uses prisms right now...
  template <class ct>
  class PrismQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    //! Constructor
    //! \param order The desired order (provided by the user)
    //! \param id A unique id (provided by QuadratureProvider)
    PrismQuadrature(const GeometryType&, int order, size_t id);

    //! The geometry type is... prism!
    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::prism,3);
    }

    //! Returns the actual order.
    virtual int order() const {
      return order_;
    }

    //! The maximal order of prism quadratures.
    static size_t maxOrder() { return PrismPoints::highest_order; }

  private:
    int order_;
  };

  //! A quadrature class for pyramids
  //! The HD stuff is used here, but needs some rework since only one rule
  //! is provided. But since nobody here uses pyramids right now...
  template <class ct>
  class PyramidQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    //! Constructor
    //! \param order The desired order (provided by the user)
    //! \param id A unique id (provided by QuadratureProvider)
    PyramidQuadrature(const GeometryType&, int order, size_t id);

    //! The geometry type is... pyramid!
    virtual GeometryType geometry() const {
      return GeometryType(GeometryType::pyramid,3);
    }

    //! Returns the actual order of the quadrature.
    virtual int order() const {
      return order_;
    }

    //! The maximal order of the pyramid quadratures.
    static size_t maxOrder() { return PyramidPoints::highest_order; }

  private:
    int order_;
  };
#endif // end USE_DUNE_QUADRATURES not defined

  //! \brief Allows injection of arbitrary points as quadrature points.
  //! Useful to test some features of the quadrature framework in isolation
  //! and with known input data. Each TestQuadrature object gets its own
  //! unique id.
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

  //! default defines for used quadratures 
  template <typename ct, int dim> struct DefaultQuadratureTraits
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory<ct,dim> CubeQuadratureType;
#else 
    typedef CubeQuadrature<ct, dim>   CubeQuadratureType; 
#endif
    typedef QuadratureImp<ct,dim>     IntegrationPointListType;
  }; 

  //! quadratures for points 
  template <typename ct> 
  struct DefaultQuadratureTraits<ct,0>  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory<ct,0> PointQuadratureType;
#else 
    typedef CubeQuadrature<ct, 0>   PointQuadratureType;     
#endif
    typedef QuadratureImp<ct,0>     IntegrationPointListType;
  };
  
  //! quadratures for lines 
  template <typename ct>
  struct DefaultQuadratureTraits<ct,1>  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory<ct,1> LineQuadratureType;
#else 
    typedef CubeQuadrature<ct, 1>   LineQuadratureType;     
#endif
    typedef QuadratureImp<ct,1>     IntegrationPointListType;
  };
  
  //! quadratures for simplex and cubes 
  template <typename ct>
  struct DefaultQuadratureTraits<ct,2>  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory<ct,2> SimplexQuadratureType;
    typedef QuadratureRulesFactory<ct,2> CubeQuadratureType;
#else
    typedef CubeQuadrature<ct, 2>    CubeQuadratureType;     
    typedef SimplexQuadrature<ct, 2> SimplexQuadratureType;     
#endif
    typedef QuadratureImp<ct,2>      IntegrationPointListType;
  };
  
  //! quadratures for simplex, cubes, prisms, and pyramids
  template <typename ct>
  struct DefaultQuadratureTraits<ct,3>  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory<ct,3> SimplexQuadratureType;
    typedef QuadratureRulesFactory<ct,3> CubeQuadratureType;

    typedef QuadratureRulesFactory<ct,3> PrismQuadratureType;
    typedef QuadratureRulesFactory<ct,3> PyramidQuadratureType;
#else 
    typedef CubeQuadrature<ct, 3>    CubeQuadratureType;     
    typedef SimplexQuadrature<ct, 3> SimplexQuadratureType;     

    typedef PrismQuadrature<ct>      PrismQuadratureType;
    typedef PyramidQuadrature<ct>    PyramidQuadratureType;
#endif

    typedef QuadratureImp<ct,3>      IntegrationPointListType;
  };
  

  //! The actual interface class for quadratures.
  //! Quadrature is a proxy for the actual implementations of the quadratures.
  //! During construction, the actual Quadrature object is configured with
  //! an appropriate implementation object from the QuadratureProvider object
  //! (Monostate pattern). The design goal here is to minimize the construction
  //! time (the actual implementations can be created once and reused as often
  //! as you like) and to insulate the user from all this initialisation and
  //! storage stuff.
  template <typename ct, int dim, 
            template <class, int> class IntegrationTraits>
  class IntegrationPointList  
  {
    typedef IntegrationPointList<ct,dim,IntegrationTraits> ThisType;
    typedef IntegrationTraits<ct,dim> Traits;
  public:
    enum { dimension = dim };

    //! type of integration point list 
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;
    //! type of coordinate 
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  public:
    //! return reference to implementation of integration point List 
    const IntegrationPointListType& ipList() const { return ipList_; }

    //! Constructor
    //! \param geo The geometry type the quadrature points belong to.
    //! \param order The order of the quadrature (i.e. polynoms up to order
    //! are integrated exactly).
    IntegrationPointList(const GeometryType& geo, int order) :
      ipList_(QuadratureProvider<ct, dim, IntegrationTraits>::getQuadrature(geo, order))
    {}

    //! Constructor for testing purposes
    //! \param quadImp Quadrature implementation for this test
    IntegrationPointList(const IntegrationPointListType& ipList) :
      ipList_(ipList)
    {}

    //! Copy constructor 
    IntegrationPointList(const IntegrationPointList& org) :
      ipList_(org.ipList_)
    {}

    //! The total number of quadrature points.
    int nop() const {
      return ipList_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(size_t i) const {
      return ipList_.point(i);
    }

    //! A unique id per quadrature type.
    //! Quadratures are considered as distinct when they differ in the
    //! following points: geometry type, order, dimension and implementation.
    //! \note At the time of writing this, there is only one implementation
    //! per geometry type, order and dimension provided, but the concept is
    //! easily extendible beyond that.
    size_t id() const {
      return ipList_.id();
    }

    //! The actual order of the quadrature.
    //! The actual order can be higher as the desired order when no 
    //! implementation for the desired order is found.
    int order() const {
      return ipList_.order();
    }

    //! The geometry type the quadrature points belong to.
    GeometryType geometry() const {
      return ipList_.geometry();
    }

  protected:
    const IntegrationPointListType& ipList_;
  };

  //! The actual interface class for quadratures.
  //! Quadrature is a proxy for the actual implementations of the quadratures.
  //! During construction, the actual Quadrature object is configured with
  //! an appropriate implementation object from the QuadratureProvider object
  //! (Monostate pattern). The design goal here is to minimize the construction
  //! time (the actual implementations can be created once and reused as often
  //! as you like) and to insulate the user from all this initialisation and
  //! storage stuff.
  template <typename ct, int dim, 
            template <class, int> class QuadratureTraits = DefaultQuadratureTraits >
  class Quadrature : public IntegrationPointList<ct,dim,QuadratureTraits> 
  {
    typedef DefaultQuadratureTraits<ct,dim> Traits; 
    typedef IntegrationPointList<ct,dim,QuadratureTraits> BaseType;
  public:
    enum { dimension = dim };

    typedef typename Traits :: IntegrationPointListType
      IntegrationPointListType;

    //! type of global coordinate vectors 
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  public:
    //! Constructor
    //! \param geo The geometry type the quadrature points belong to.
    //! \param order The order of the quadrature (i.e. polynoms up to order
    //! are integrated exactly).
    Quadrature(const GeometryType& geo, int order) :
      BaseType(geo,order)
    {}

    //! Constructor for testing purposes
    //! \param quadImp Quadrature implementation for this test
    Quadrature(const IntegrationPointListType& ipList) :
      BaseType(ipList)
    {}

    //! Copy constructor
    Quadrature(const Quadrature& org) :
      BaseType(org)
    {}

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const ct& weight(size_t i) const {
      return this->ipList().weight(i);
    }
  };

} // end namespace Dune

#include "quadrature.cc"

#undef USE_DUNE_QUADRATURES
#undef HAVE_ALBERTA_FOUND

#endif
