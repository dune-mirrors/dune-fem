#ifndef DUNE_ADIQUADRATURE_HH
#define DUNE_ADIQUADRATURE_HH

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

// undef HAVE_ALBERTA_FOUND, because quaratures boogie  
#undef HAVE_ALBERTA_FOUND

#ifdef HAVE_ALBERTA_FOUND
// inlcude albertagrid.hh includes the needed alberta.h 
#include <dune/grid/albertagrid.hh>
#endif

#include "quadprovider.hh"
#include "simplexpoints.hh"
#include "gausspoints.hh"
#include "prismpoints.hh"
#include "pyramidpoints.hh"
#include "dunequadratures.hh"

namespace Dune {

  // Forward declaration
  template <typename ct, int dim, template <class,int> class QuadratureTraits>  
  class QuadratureProvider;
  template <class ct, int dim>
  class QuadratureImp;
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

  //! Generic implementation of a quadrature.
  //! A quadrature in the Dune sense is nothing but a set of points and
  //! weights.
  template <typename ct, int dim>
  class IntegrationPointListImp {
  public:
    //! Local coordinate type for the quadrature.
    typedef FieldVector<ct, dim> CoordinateType;
    
    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  public:
    //! Virtual destructor
    virtual ~IntegrationPointListImp() {}

    //! Coordinates of integration point i.
    const CoordinateType& point(size_t i) const {
      return points_[i];
    }

    //! Number of integration points
    size_t nop() const {
      return points_.size();
    }

    //! A globally unique identifier amongst all integration point lists 
    //! (even for other dimensions).
    size_t id() const {
      return id_;
    }

    //! Maximal degree of polynomial that gets integrated exactly by the
    //! quadrature.
    virtual int order() const = 0;

    //! Geometry type the integration point list is defined for.
    virtual GeometryType geometry() const = 0;

  protected:
    //! QuadratureImps are filled by the derived classes
    IntegrationPointListImp(size_t id);

    //! Adds a quadrature point/weight pair. To be used in the constructor
    //! of derived class.
    void addIntegrationPoint(const CoordinateType& point);

  private:
    //! Copying is forbidden!
    IntegrationPointListImp(const IntegrationPointListImp&);

  protected:
    std::vector<CoordinateType> points_;
    const size_t id_;
  };

  template <typename ct, int dim>
  class QuadratureImp : public IntegrationPointListImp<ct,dim> {
    // type of base class 
    typedef IntegrationPointListImp<ct,dim> BaseType;
  public:
    //! Local coordinate type for the quadrature.
    typedef typename BaseType :: CoordinateType  CoordinateType;
    
  public:
    //! Virtual destructor
    virtual ~QuadratureImp() {}

    //! The weight of quadrature point i.
    const ct& weight(size_t i) const {
      return weights_[i];
    }

  protected:
    //! QuadratureImps are filled by the derived classes
    QuadratureImp(size_t id);

    //! Adds a quadrature point/weight pair. To be used in the constructor
    //! of derived class.
    void addQuadraturePoint(const CoordinateType& point, ct weight);

  private:
    //! Copying is forbidden!
    QuadratureImp(const QuadratureImp&);

  private:
    // vector holding weights of each integration point 
    std::vector<ct> weights_;
  };

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
    typedef CubeQuadrature<ct, dim> CubeQuadratureType; 
    typedef QuadratureImp<ct,dim>     IntegrationPointListType;
  }; 

  //! quadratures for points 
  template <typename ct> 
  struct DefaultQuadratureTraits<ct,0>  
  {
    typedef CubeQuadrature<ct, 0>   PointQuadratureType;     
    typedef QuadratureImp<ct,0>     IntegrationPointListType;
  };
  
  //! quadratures for lines 
  template <typename ct>
  struct DefaultQuadratureTraits<ct,1>  
  {
    typedef CubeQuadrature<ct, 1>   LineQuadratureType;     
    typedef QuadratureImp<ct,1>     IntegrationPointListType;
  };
  
  //! quadratures for simplex and cubes 
  template <typename ct>
  struct DefaultQuadratureTraits<ct,2>  
  {
    typedef CubeQuadrature<ct, 2>    CubeQuadratureType;     
    typedef SimplexQuadrature<ct, 2> SimplexQuadratureType;     
    /*
    typedef QuadratureRulesFactory<ct,2> SimplexQuadratureType;
    typedef QuadratureRulesFactory<ct,2> CubeQuadratureType;
    */
    typedef QuadratureImp<ct,2>      IntegrationPointListType;
  };
  
  //! quadratures for simplex, cubes, prisms, and pyramids
  template <typename ct>
  struct DefaultQuadratureTraits<ct,3>  
  {
    typedef CubeQuadrature<ct, 3>    CubeQuadratureType;     
    typedef SimplexQuadrature<ct, 3> SimplexQuadratureType;     

    typedef PrismQuadrature<ct>      PrismQuadratureType;
    typedef PyramidQuadrature<ct>    PyramidQuadratureType;

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
#undef HAVE_ALBERTA_FOUND

#endif
