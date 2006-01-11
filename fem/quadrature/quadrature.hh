#ifndef DUNE_ADIQUADRATURE_HH
#define DUNE_ADIQUADRATURE_HH

#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#if HAVE_ALBERTA 
// inlcude albertagrid.hh includes the needed alberta.h 
#include <dune/grid/albertagrid.hh>
#endif

#include "quadprovider.hh"
#include "simplexpoints.hh"
#include "gausspoints.hh"
#include "prismpoints.hh"
#include "pyramidpoints.hh"

namespace Dune {

  // Forward declaration
  template <typename ct, int dim>
  class QuadratureProvider;

  //! Generic implementation of a quadrature.
  //! A quadrature in the Dune sense is nothing but a set of points and
  //! weights.
  template <typename ct, int dim>
  class QuadratureImp {
  public:
    //! Local coordinate type for the quadrature.
    typedef FieldVector<ct, dim> CoordinateType;
    
  public:
    //! Virtual destructor
    virtual ~QuadratureImp() {}

    //! Coordinates of quadrature point i.
    const CoordinateType& point(size_t i) const {
      return points_[i];
    }

    //! The weight of quadrature point i.
    const ct& weight(size_t i) const {
      return weights_[i];
    }

    //! Number of quadrature points
    size_t nop() const {
      return points_.size();
    }

    //! A globally unique identifier amongst all quadratures (even for other
    //! dimensions).
    size_t id() const {
      return id_;
    }

    //! Maximal degree of polynomial that gets integrated exactly by the
    //! quadrature.
    virtual int order() const = 0;

    //! Geometry type the quadrature defines points for.
    virtual GeometryType geometry() const = 0;

  protected:
    //! QuadratureImps are filled by the derived classes
    QuadratureImp(size_t id);

    //! Adds a quadrature point/weight pair. To be used in the constructor
    //! of derived class.
    void addQuadraturePoint(const CoordinateType& point, ct weight);

  private:
    //! Copying is forbidden!
    QuadratureImp(const QuadratureImp&);
    //! Assignment is forbidden!
    //QuadratureImp operator=(const QuadratureImp&);

  private:
    std::vector<CoordinateType> points_;
    std::vector<ct> weights_;

    size_t id_;
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
#if HAVE_ALBERTA 
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
    SimplexQuadrature(int order, size_t id);
    
    //! The geometry type is... simplex!
    virtual GeometryType geometry() const {
      return simplex;
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
    CubeQuadrature(int order, size_t id);
    
    //! The geometry type is... cube!
    virtual GeometryType geometry() const {
      return cube;
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
    LineQuadrature(int order, size_t id);

    virtual GeometryType geometry() const {
      return cube;
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
    TriangleQuadrature(int order, size_t id);

    virtual GeometryType geometry() const {
      return simplex;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { 
#if HAVE_ALBERTA 
      return 7; 
#else 
      return 5; 
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
    QuadrilateralQuadrature(int order, size_t id);

    virtual GeometryType geometry() const {
      return cube;
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
    TetraQuadrature(int order, size_t id);

    virtual GeometryType geometry() const {
      return simplex;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { 
#if HAVE_ALBERTA 
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

  //! A quadrature class for hexahedra
  //! \note This class is redundant as CubeQuadrature can be used instead.
  template <class ct>
  class HexaQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    HexaQuadrature(int order, size_t id);

    virtual GeometryType geometry() const {
      return cube;
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
    PrismQuadrature(int order, size_t id);

    //! The geometry type is... prism!
    virtual GeometryType geometry() const {
      return prism;
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
    PyramidQuadrature(int order, size_t id);

    //! The geometry type is... pyramid!
    virtual GeometryType geometry() const {
      return pyramid;
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
    TestQuadrature(GeometryType geo, int order);

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


  //! The actual interface class for quadratures.
  //! Quadrature is a proxy for the actual implementations of the quadratures.
  //! During construction, the actual Quadrature object is configured with
  //! an appropriate implementation object from the QuadratureProvider object
  //! (Monostate pattern). The design goal here is to minimize the construction
  //! time (the actual implementations can be created once and reused as often
  //! as you like) and to insulate the user from all this initialisation and
  //! storage stuff.
  template <typename ct, int dim>
  class Quadrature 
  {
  public:
    enum { dimension = dim };

    typedef typename QuadratureImp<ct, dim>::CoordinateType CoordinateType;

  public:
    //! Constructor
    //! \param geo The geometry type the quadrature points belong to.
    //! \param order The order of the quadrature (i.e. polynoms up to order
    //! are integrated exactly).
    Quadrature(GeometryType geo, int order) :
      quad_(QuadratureProvider<ct, dim>::getQuadrature(geo, order))
    {}

    //! Constructor for testing purposes
    //! \param quadImp Quadrature implementation for this test
    Quadrature(const QuadratureImp<ct, dim>& quadImp) :
      quad_(quadImp)
    {}

    //! The total number of quadrature points.
    int nop() const {
      return quad_.nop();
    }

    //! Access to the ith quadrature point.
    const CoordinateType& point(size_t i) const {
      return quad_.point(i);
    }

    //! Access to the weight of quadrature point i.
    //! The quadrature weights sum up to the volume of the respective reference
    //! element.
    const ct& weight(size_t i) const {
      return quad_.weight(i);
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
      return quad_.geometry();
    }

  private:
    const QuadratureImp<ct, dim>& quad_;
  };

} // end namespace Dune

#include "include_cc.hh"
#include "quadrature.cc"

#endif
