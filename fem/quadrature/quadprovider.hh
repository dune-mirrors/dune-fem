#ifndef DUNE_QUADPROVIDER_HH
#define DUNE_QUADPROVIDER_HH

#include "quadrature.hh"
#include "idprovider.hh"

namespace Dune {

  // Forward declarations
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

  template <typename ct, int dim>
  class QuadratureProvider {
  public:
    static const QuadratureImp<ct, dim>& getQuadrature(GeometryType geo, 
                                                       int order);
  };

  // Specialisaions
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
    static std::vector<SimplexQuadrature<ct, 2>*> triangleQuads_;
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
    static std::vector<SimplexQuadrature<ct, 3>*> tetraQuads_;
    static std::vector<CubeQuadrature<ct, 3>*> hexaQuads_;
    static std::vector<PrismQuadrature<ct>*> prismQuads_;
    static std::vector<PyramidQuadrature<ct>*> pyramidQuads_;
  };

#include "quadprovider.cc"
}

#endif
