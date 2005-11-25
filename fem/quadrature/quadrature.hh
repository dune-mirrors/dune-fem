#ifndef DUNE_ADIQUADRATURE_HH
#define DUNE_ADIQUADRATURE_HH

#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include "quadprovider.hh"
#include "simplexpoints.hh"
#include "gausspoints.hh"
#include "prismpoints.hh"
#include "pyramidpoints.hh"

namespace Dune {

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
    virtual GeometryType geo() const = 0;

  protected:
    //! QuadratureImps are filled by the derived classes
    QuadratureImp(size_t id);

    //! Adds a quadrature point/weight pair. To be used in the constructor
    //! of derived class.
    void addQuadraturePoint(const CoordinateType& point, ct weight);

  private:
    QuadratureImp(const QuadratureImp&);
    QuadratureImp operator=(const QuadratureImp&);

  private:
    std::vector<CoordinateType> points_;
    std::vector<ct> weights_;

    size_t id_;
  };

  template <class ct, int dim>
  class SimplexQuadrature : public QuadratureImp<ct, dim>
  {
  public:
    typedef FieldVector<ct, dim> CoordinateType;

  public:
    SimplexQuadrature(int order, size_t id);
    
    virtual GeometryType geo() const {
      return cube;
    }
    
    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return 19; }

  private:
    int order_;
  };

  template <class ct, int dim>
  class CubeQuadrature : public QuadratureImp<ct, dim>
  {
  public:
    typedef FieldVector<ct, dim> CoordinateType;

  public:
    CubeQuadrature(int order, size_t id);
    
    virtual GeometryType geo() const {
      return cube;
    }
    
    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPoints::highestOrder; }

  private:
    int order_;
  };

  template <class ct>
  class LineQuadrature : public QuadratureImp<ct, 1> 
  {
  public:
    typedef FieldVector<ct, 1> CoordinateType;

  public:
    LineQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return cube;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPoints::highestOrder; }

  private:
    int order_;
  };

  template <class ct>
  class TriangleQuadrature : public QuadratureImp<ct, 2>
  {
  public:
    typedef FieldVector<ct, 2> CoordinateType;

  public:
    TriangleQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return simplex;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return 0; }


  private:
    int order_;
  };

  template <class ct>
  class QuadrilateralQuadrature : public QuadratureImp<ct, 2>
  {
  public:
    typedef FieldVector<ct, 2> CoordinateType;

  public:
    QuadrilateralQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return cube;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPoints::highestOrder; }

  private:
    int order_;
  };

  template <class ct>
  class TetraQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    TetraQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return simplex;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return 0; }

  private:
    int order_;
  };

  template <class ct>
  class HexaQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    HexaQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return cube;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return GaussPoints::highestOrder; }

  private:
    int order_;
  };

  template <class ct>
  class PrismQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    PrismQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return prism;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return PrismPoints::highest_order; }

  private:
    int order_;
  };

  template <class ct>
  class PyramidQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    PyramidQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return pyramid;
    }

    virtual int order() const {
      return order_;
    }

    static size_t maxOrder() { return PyramidPoints::highest_order; }

  private:
    int order_;
  };



  template <typename ct, int dim>
  class Quadrature 
  {
  public:
    enum { dimension = dim };

    typedef typename QuadratureImp<ct, dim>::CoordinateType CoordinateType;

  public:
    Quadrature(GeometryType geo, int order) :
      quad_(QuadratureProvider<ct, dim>::getQuadrature(geo, order))
    {}

    const CoordinateType& point(size_t i) const {
      return quad_.point(i);
    }

    int nop() const {
      return quad_.nop();
    }

    const ct& weight(size_t i) const {
      return quad_.weight(i);
    }

    size_t id() const {
      return quad_.id();
    }

    int order() const {
      return quad_.order();
    }

    GeometryType geo() const {
      return quad_.geo();
    }
  private:
    const QuadratureImp<ct, dim>& quad_;
  };

} // end namespace Dune

#include "quadrature.cc"

#endif
