#include <config.h>

#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/utility/twistutility.hh>

#include "cachequad_test.hh"
#include "albertagrid_fixture.hh"
#include "alugrid_fixture.hh"
#include "sgrid_fixture.hh"

namespace Dune {

  void CachingQuadrature_Test::run() 
  {
    codim0Test();
    codim1AlbertaTest();
    codim1ALUHexaTest();
    codim1ALUTetraTest();
    codim1SGridTest();
  }

  void CachingQuadrature_Test::codim0Test() 
  {
    const int dim = 2;
    const int codim = 0;

    typedef AlbertaGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;

    GridFixtureType fix(albertaGridFile_);
    GridType& grid = fix.grid();
    
    QuadratureType quad(*grid.leafbegin<0>(), 3);
    const PointVectorType& points = 
      PointProviderType::getPoints(quad.id(), triangle);

    _test(points.size() == (size_t)quad.nop());
    
    for (int i = 0; i < quad.nop(); ++i) {
      for (int d = 0; d < dim; ++d) {
        _floatTest(points[i][d], quad.point(quad.cachingPoint(i))[d]);
      }
    }
  }

  void CachingQuadrature_Test::codim1AlbertaTest() 
  {
    const int dim = 2;
    const int codim = 1;

    typedef AlbertaGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;

    GridFixtureType fix(albertaGridFile_);
    GridType& grid = fix.grid();

    IntersectionIterator endit = grid.leafbegin<0>()->iend();
    IntersectionIterator it = grid.leafbegin<0>()->ibegin();
    for (; it != endit; ++it) {
      const LocalGeometryType& geo = it.intersectionSelfLocal();
      QuadratureType quad(it, 4, 0, ElementQuadrature<GridType, 1>::INSIDE);
      const PointVectorType& points = 
        PointProviderType::getPoints(quad.id(), triangle);

      _test(points.size() == (size_t)3*quad.nop());
      
      for (int i = 0; i < quad.nop(); ++i) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
        }
      }
    }
  }

  
  void CachingQuadrature_Test::codim1ALUHexaTest() 
  {
    const int dim = 3;
    const int codim = 1;

    typedef ALUGridFixture<hexa> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;

    GridFixtureType fix(aluGridHexaFile_);
    GridType& grid = fix.grid();

    IntersectionIterator endit = grid.leafbegin<0>()->iend();
    IntersectionIterator it = grid.leafbegin<0>()->ibegin();
    for (; it != endit; ++it) {
      const LocalGeometryType& geo = it.intersectionSelfLocal();
      TwistUtility<GridType> twistUtil(grid);
      int twist = twistUtil.twistInSelf(it);
      //std::cout << "Twist == " << twist << std::endl;
      QuadratureType quad(it, 4, twist, ElementQuadrature<GridType, 1>::INSIDE);
      const PointVectorType& points = 
        PointProviderType::getPoints(quad.id(), cube);

      _test(points.size() == (size_t)6*quad.nop());
      
      for (int i = 0; i < quad.nop(); ++i) {
        for (int d = 0; d < dim; ++d) {
          // geo has twist
          _floatTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
        }
        //std::cout << "pt " << i << ": " << points[quad.cachingPoint(i)]
        //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
      }
      
    } // end iterator loop
  }

  void CachingQuadrature_Test::codim1ALUTetraTest() 
  {
    const int dim = 3;
    const int codim = 1;

    typedef ALUGridFixture<tetra> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;

    GridFixtureType fix(aluGridTetraFile_);
    GridType& grid = fix.grid();

    IntersectionIterator endit = grid.leafbegin<0>()->iend();
    IntersectionIterator it = grid.leafbegin<0>()->ibegin();
    for (; it != endit; ++it) {
      const LocalGeometryType& geo = it.intersectionSelfLocal();
      TwistUtility<GridType> twistUtil(grid);
      int twist = twistUtil.twistInSelf(it);
      //std::cout << "Twist == " << twist << std::endl;
      QuadratureType quad(it, 4, twist, ElementQuadrature<GridType, 1>::INSIDE);
      const PointVectorType& points = 
        PointProviderType::getPoints(quad.id(), simplex);

      _test(points.size() == (size_t)4*quad.nop());
      
      for (int i = 0; i < quad.nop(); ++i) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
        }
        //std::cout << "pt " << i << ": " << points[quad.cachingPoint(i)]
        //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
      }
      
    } // end iterator loop
  }

  void CachingQuadrature_Test::codim1SGridTest() 
  {
    const int dim = 2;
    const int codim = 1;

    typedef SGridFixture<dim, dim> GridFixtureType;
    typedef GridFixtureType::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef CachingQuadrature<GridType, codim> QuadratureType;
    typedef PointProvider<GridType::ctype, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef IntersectionIterator::LocalGeometry LocalGeometryType;

    const int N = 2;

    GridFixtureType fix(N);
    GridType& grid = fix.grid();

  _fail("if you set the order to 4 here, the cache provider will find the wrong quadrature, since there is no distinction between line quadrature for triangles and quadrilaterals right now");


    IntersectionIterator endit = grid.leafbegin<0>()->iend();
    IntersectionIterator it = grid.leafbegin<0>()->ibegin();
    for (; it != endit; ++it) {
      const LocalGeometryType& geo = it.intersectionSelfLocal();
      TwistUtility<GridType> twistUtil(grid);
      int twist = twistUtil.twistInSelf(it);
    
      _test(twist == 0);

      QuadratureType quad(it, 5, twist, ElementQuadrature<GridType, 1>::INSIDE);
      const PointVectorType& points = 
        PointProviderType::getPoints(quad.id(), cube);

      _test(points.size() == (size_t)2*dim*quad.nop());
      
      for (int i = 0; i < quad.nop(); ++i) {
        for (int d = 0; d < dim; ++d) {
          _floatTest(points[quad.cachingPoint(i)][d],
                     geo.global(quad.localPoint(i))[d]);
        }
        //std::cout << "pt " << i << ": " << points[quad.cachingPoint(i)]
        //          << " == " << geo.global(quad.localPoint(i)) << std::endl;
      }
      
    } // end iterator loop
  }


} // end namespace Dune
 
