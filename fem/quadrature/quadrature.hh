#ifndef DUNE_QUADRATURE_HH
#define DUNE_QUADRATURE_HH

#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include "idprovider.hh"

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

  private:
    int order_;
  };

  template <class ct>
  class TetrahedronQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    TetrahedronQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return simplex;
    }

    virtual int order() const {
      return order_;
    }

  private:
    int order_;
  };

  template <class ct>
  class HexahedronQuadrature : public QuadratureImp<ct, 3>
  {
  public:
    typedef FieldVector<ct, 3> CoordinateType;

  public:
    HexahedronQuadrature(int order, size_t id);

    virtual GeometryType geo() const {
      return cube;
    }

    virtual int order() const {
      return order_;
    }

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

  // Specialisations
  template <typename ct>
  class QuadratureProvider<ct, 1>
  {
  public:
    static const QuadratureImp<ct, 1>& getQuadrature(GeometryType geo, 
                                                     int order) {
      assert(geo == cube  || geo == simplex || geo == line);
      assert(order >= 0);

      if (!quads_[order]) {
        makeQuad(order);
      }

      return *quads_[order];
    }
  private:
    QuadratureProvider();
    QuadratureProvider(const QuadratureProvider&);
    QuadratureProvider& operator=(const QuadratureProvider&);

    static void makeQuad(int order) {
      quads_.resize(std::max(order+1, quads_.size()));
      quads_[order] = 
        new LineQuadrature<ct>(order, IdProvider::instance().newId());
    }
    
  private:
    static std::vector<QuadratureImp<ct, 1>*> quads_;
  }; 

  template <typename ct>
  std::vector<QuadratureImp<ct, 1>*> QuadratureProvider<ct, 1>::quads_;

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
        if (!triangleQuads_[order]) {
          makeTriangleQuad(order);
        }
        return triangleQuads_[order];
      case cube:
      case quadrilateral:
        if (!quadrilateralQuads_[order]) {
          makeQuadrilateralQuad(order);
        }
        return quadrilateralQuads_[order];
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

    static void makeTriangleQuad(int order) {
      triangleQuads_.resize(std::max(order+1, triangleQuads_.size()));
      triangleQuads_[order] = 
        new TriangleQuadrature<ct>(order, IdProvider::instance().newId());
    }

    static void makeQuadrilateralQuad(int order) {
      quadrilateralQuads_.resize(std::max(order+1,quadrilateralQuads_.size()));
      quadrilateralQuads_[order] = 
        new QuadrilateralQuadrature<ct>(order, IdProvider::instance().newId());
    }

  private:
    static std::vector<QuadratureImp<ct, 2>*> triangleQuads_;
    static std::vector<QuadratureImp<ct, 2>*> quadrilateralQuads_;
  };

  template <typename ct>
  std::vector<QuadratureImp<ct, 2>*> QuadratureProvider<ct, 2>::triangleQuads_;
  template <typename ct>
  std::vector<QuadratureImp<ct, 2>*> QuadratureProvider<ct, 2>::quadrilateralQuads_;

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
        if (!tetraQuads_[order]) {
          makeTetrahedronQuad(order);
        }
        return tetraQuads_[order];
      case cube:
      case hexahedron:
        if (!hexaQuads_[order]) {
          makeHexahedronQuad(order);
        }
        return hexaQuads_[order];
      case prism:
        if (!prismQuads_[order]) {
          makePrismQuad(order);
        }
        return prismQuads_[order];
      case pyramid:
        if (!pyramidQuads_[order]) {
          makePyramidQuad(order);
        }
        return pyramidQuads_[order];
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

    static void makeTetrahedronQuad(int order) {
      tetraQuads_.resize(std::max(order+1, tetraQuads_.size()));
      tetraQuads_[order] = 
        new TetrahedronQuadrature<ct>(order, IdProvider::instance().newId());
    }

    static void makeHexahedronQuad(int order) {
      hexaQuads_.resize(std::max(order+1, hexaQuads_.size()));
      hexaQuads_[order] = 
        new HexahedronQuadrature<ct>(order, IdProvider::instance().newId());
    }

    static void makePrismQuad(int order) {
      prismQuads_.resize(std::max(order+1, prismQuads_.size()));
      prismQuads_[order] = 
        new PrismQuadrature<ct>(order, IdProvider::instance().newId());
    }

    static void makePyramidQuad(int order) {
      pyramidQuads_.resize(std::max(order+1, pyramidQuads_.size()));
      pyramidQuads_[order] = 
        new PyramidQuadrature<ct>(order, IdProvider::instance().newId());
    }

  private:
    static std::vector<QuadratureImp<ct, 3>*> tetraQuads_;
    static std::vector<QuadratureImp<ct, 3>*> hexaQuads_;
    static std::vector<QuadratureImp<ct, 3>*> prismQuads_;
    static std::vector<QuadratureImp<ct, 3>*> pyramidQuads_;
  };

  template <typename ct>
  std::vector<QuadratureImp<ct, 3>*> QuadratureProvider<ct, 3>::tetraQuads_;
  template <typename ct>
  std::vector<QuadratureImp<ct, 3>*> QuadratureProvider<ct, 3>::hexaQuads_;
  template <typename ct>
  std::vector<QuadratureImp<ct, 3>*> QuadratureProvider<ct, 3>::prismQuads_;
  template <typename ct>
  std::vector<QuadratureImp<ct, 3>*> QuadratureProvider<ct, 3>::pyramidQuads_;
} // end namespace Dune

#include "quadrature.cc"

#endif
