#include <config.h>

#include <dune/grid/sgrid.hh>
#include "sgrid_fixture.hh"

#ifdef ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>
#include "alugrid_fixture.hh"
#endif

#ifdef ENABLE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include "albertagrid_fixture.hh"
#endif

#include "cache_test.hh"

namespace Dune {
  
  void CacheProvider_Test::run() 
  {
    hexaTest();
    tetraTest();
    triangleTest();
    quadTest();
  }

  void CacheProvider_Test::hexaTest()
  {
#ifdef ENABLE_ALUGRID
    const int dim = 3;
    const int codim = 1;

    typedef ALUGridFixture<hexa> FixType;
    typedef FixType::GridType GridType;
    typedef CacheProvider<GridType, codim> CacheProviderType;
    typedef Quadrature< GridType :: ctype, GridType :: dimension-1> QuadratureType;
    //typedef CacheProviderType::QuadratureType QuadratureType;
    typedef CacheProviderType::MapperType MapperType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GeometryType elemGeo = GeometryType(GeometryType::cube,3);
    GeometryType faceGeo = GeometryType(GeometryType::cube,2);

    // Get reference element
    const ReferenceElement<double, dim>& refElem =
      ReferenceElements<double, dim>::general(elemGeo);
    
    // Build quadrature
    QuadratureType quad(faceGeo, 3);

    // Ask for one mapper so that the points get registered
    CacheProviderType::getMapper(quad, elemGeo, 0, 0);

    const PointVectorType& points =
      PointProviderType::getPoints(quad.id(), elemGeo);

    // Loop over all faces
    for (int i = 0; i < refElem.size(codim); ++i) {
      const MapperType& m = 
        CacheProviderType::getMapper(quad, elemGeo, i, 0);
      
      _test(m.size() == (size_t) quad.nop());
      // Loop over all points
      for (size_t j = 0; j < m.size(); ++j) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[m[j]][d], 
                     refElem.global<codim>(quad.point(j), i, codim)[d]);
        }
      }

    }
#endif
  }

  void CacheProvider_Test::tetraTest()
  {
#ifdef ENABLE_ALUGRID
    const int dim = 3;
    const int codim = 1;

    typedef ALUGridFixture<tetra> FixType;
    typedef FixType::GridType GridType;
    typedef CacheProvider<GridType, codim> CacheProviderType;
    //typedef CacheProviderType::QuadratureType QuadratureType;
    typedef Quadrature< GridType :: ctype, GridType :: dimension-1> QuadratureType;
    typedef CacheProviderType::MapperType MapperType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GeometryType elemGeo = GeometryType(GeometryType::simplex,3);
    GeometryType faceGeo = GeometryType(GeometryType::simplex,2);

    // Get reference element
    const ReferenceElement<double, dim>& refElem =
      ReferenceElements<double, dim>::general(elemGeo);
    
    // Build quadrature
    QuadratureType quad(faceGeo, 3);

    // Ask for one mapper so that the points get registered
    CacheProviderType::getMapper(quad, elemGeo, 0, 0);

    const PointVectorType& points =
      PointProviderType::getPoints(quad.id(), elemGeo);

    // Loop over all faces
    for (int i = 0; i < refElem.size(codim); ++i) {
      const MapperType& m = CacheProviderType::getMapper(quad, elemGeo, i, 0);
      
      _test(m.size() == (size_t) quad.nop());
      // Loop over all points
      for (size_t j = 0; j < m.size(); ++j) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[m[j]][d], 
                     refElem.global<codim>(quad.point(j), i, codim)[d]);
        }
      }
    }
#endif
  }

  void CacheProvider_Test::triangleTest()
  {
#ifdef ENALBE_ALBERTA
    const int dim = 2;
    const int codim = 1;

    typedef AlbertaGridFixture<dim, dim> FixType;
    typedef FixType::GridType GridType;
    typedef CacheProvider<GridType, codim> CacheProviderType;
    //typedef CacheProviderType::QuadratureType QuadratureType;
    typedef Quadrature< GridType :: ctype, GridType :: dimension-1> QuadratureType;
    typedef CacheProviderType::MapperType MapperType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GeometryType elemGeo = GeometryType(GeometryType::simplex,2);
    GeometryType faceGeo = GeometryType(GeometryType::simplex,1);

    // Get reference element
    const ReferenceElement<double, dim>& refElem =
      ReferenceElements<double, dim>::general(elemGeo);
    
    // Build quadrature
    QuadratureType quad(faceGeo, 3);

    // Ask for one mapper so that the points get registered
    CacheProviderType::getMapper(quad, elemGeo, 0, 0);

    const PointVectorType& points =
      PointProviderType::getPoints(quad.id(), elemGeo);

    // Loop over all faces
    for (int i = 0; i < refElem.size(codim); ++i) {
      const MapperType& m = CacheProviderType::getMapper(quad, elemGeo, i, 0);
      
      _test(m.size() == (size_t) quad.nop());
      // Loop over all points
      for (size_t j = 0; j < m.size(); ++j) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[m[j]][d], 
                     refElem.global<codim>(quad.point(j), i, codim)[d]);
        }
      }
    }
#endif
  }

  void CacheProvider_Test::quadTest()
  {
    const int dim = 2;
    const int codim = 1;

    typedef SGridFixture<dim, dim> FixType;
    typedef FixType::GridType GridType;
    typedef CacheProvider<GridType, codim> CacheProviderType;
    typedef Quadrature< GridType :: ctype, GridType :: dimension-1> QuadratureType;
    //typedef CacheProviderType::QuadratureType QuadratureType;
    typedef CacheProviderType::MapperType MapperType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GeometryType elemGeo = GeometryType(GeometryType::cube,2);
    GeometryType faceGeo = GeometryType(GeometryType::cube,1);

    // Get reference element
    const ReferenceElement<double, dim>& refElem =
      ReferenceElements<double, dim>::general(elemGeo);
    
    // Build quadrature
    QuadratureType quad(faceGeo, 5);

    // Ask for one mapper so that the points get registered
    CacheProviderType::getMapper(quad, elemGeo, 0, 0);

    const PointVectorType& points =
      PointProviderType::getPoints(quad.id(), elemGeo);

    // Loop over all faces
    for (int i = 0; i < refElem.size(codim); ++i) {
      const MapperType& m = CacheProviderType::getMapper(quad, elemGeo, i, 0);
      
      _test(m.size() == (size_t) quad.nop());
      // Loop over all points
      for (size_t j = 0; j < m.size(); ++j) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[m[j]][d], 
                     refElem.global<codim>(quad.point(j), i, codim)[d]);
        }
      }
    }
  }

} // end namespace Dune
