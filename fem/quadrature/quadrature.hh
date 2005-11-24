#ifndef DUNE_QUADRATURE_HH
#define DUNE_QUADRATURE_HH

#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include "idprovider.hh"
#include "gausspoints.hh"

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

    static int maxOrder() { return GaussPoints::highestOrder; }

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

    static int maxOrder() { return GaussPoints::highestOrder; }

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

    static int maxOrder() { return 0; }


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

    static int maxOrder() { return GaussPoints::highestOrder; }

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

    static int maxOrder() { return 0; }

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

    static int maxOrder() { return GaussPoints::highestOrder; }

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

    static int maxOrder() { return 0; }

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

    static int maxOrder() { return 0; }

  private:
    int order_;
  };



  template <typename ct, int dim>
  class Quadrature 
  {
  public:
    typedef typename QuadratureImp<ct, dim>::CoordinateType CoordinateType;

  public:
    Quadrature(GeometryType geo, int order) :
      quad_(QuadratureProvider<ct, dim>::getQuadrature(geo, order))
    {}

    const CoordinateType& point(size_t i) const {
      return quad_.point(i);
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

  template <typename ct, int dim>
  class QuadratureProvider {
  public:
    static const QuadratureImp<ct, dim>& getQuadrature(GeometryType geo, 
                                                       int order);
  };

  // Specialisaions

  class QuadCreator {
  public:
  template <class QuadImp>
    static const QuadImp& provideQuad(int order, std::vector<QuadImp*>& vec) {
      assert(vec.size() > static_cast<size_t>(order));
      if (!vec[order]) {
        vec[order] = new QuadImp(order, IdProvider::instance().newId());
      }
      return *vec[order];
    }
  };

  template <typename ct>
  class QuadratureProvider<ct, 1> 
  {
  public:
    static const QuadratureImp<ct, 1>& getQuadrature(GeometryType geo, 
                                                     int order) {
      assert(geo == cube  || geo == simplex || geo == line);
      assert(order >= 0);

      return QuadCreator::provideQuad(order, quads_);
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);

  private:
    static std::vector<CubeQuadrature<ct, 1>*> quads_;
  }; 

  template <typename ct>
  class QuadratureProvider<ct, 2> 
  {
  public:
    static const QuadratureImp<ct, 2>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert(geo == triangle || geo == quadrilateral || 
             geo == cube || geo == simplex);
      assert(order >= 0);

      switch (geo) {
      case simplex:
      case triangle:
        return QuadCreator::provideQuad(order, triangleQuads_);
      case cube:
      case quadrilateral:
        return QuadCreator::provideQuad(order, quadrilateralQuads_);
      default:
        DUNE_THROW(RangeError, "Element type not available for dim == 2");
      }
      // dummy return
      return *triangleQuads_[0]; 
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);

  private:
    static std::vector<TriangleQuadrature<ct>*> triangleQuads_;
    static std::vector<CubeQuadrature<ct, 2>*> quadrilateralQuads_;
  };

  template <class ct>
  class QuadratureProvider<ct, 3> 
  {
  public:
    static const QuadratureImp<ct, 3>& getQuadrature(GeometryType geo,
                                                     int order) {
      assert(geo == cube || geo == simplex || geo == tetrahedron ||
             geo == tetrahedron || geo == prism || geo == pyramid);
      assert(order >= 0);

      switch (geo) {
      case simplex:
      case tetrahedron:
        return QuadCreator::provideQuad(order, tetraQuads_);
      case cube:
      case hexahedron:
        return QuadCreator::provideQuad(order, hexaQuads_);
      case prism:
        return QuadCreator::provideQuad(order, prismQuads_);
      case pyramid:
        return QuadCreator::provideQuad(order, pyramidQuads_);
      default:
        DUNE_THROW(RangeError, "Element type not available for dim == 3");
      }
      // dummy return
      return *tetraQuads_[0];
    }

  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);

  private:
    static std::vector<TetraQuadrature<ct>*> tetraQuads_;
    static std::vector<CubeQuadrature<ct, 3>*> hexaQuads_;
    static std::vector<PrismQuadrature<ct>*> prismQuads_;
    static std::vector<PyramidQuadrature<ct>*> pyramidQuads_;
  };


} // end namespace Dune

#include "quadrature.cc"

#endif
