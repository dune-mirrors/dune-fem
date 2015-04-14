#include <config.h>
#include <set>

#include "quad_test.hh"

namespace Dune {
  namespace Fem {

  namespace {
    static void doTest( const double& a, const double& b )
    {
      if( std::abs( a - b ) > 1e-12 )
      {
        assert( false );
        std::abort();
      }
    }

    static void doTest( const bool value )
    {
      if( ! value )
      {
        assert( false );
        std::abort();
      }
    }
  }

  void Quad_Test::run() {

    // 1d types
    GeometryType line ( GeometryType::cube, 1);
    GeometryType cube1 = line;
    GeometryType simplex1 = line;


    // 2d types
    GeometryType quadrilateral ( GeometryType::cube, 2);
    GeometryType cube2 = quadrilateral;
    GeometryType triangle ( GeometryType::simplex, 2);
    GeometryType simplex2 = triangle;

    // 3d types
    GeometryType hexahedron ( GeometryType::cube, 3);
    GeometryType cube3 = hexahedron;
    GeometryType tetrahedron ( GeometryType::simplex, 3);
    GeometryType simplex3 = tetrahedron;
    GeometryType prism ( GeometryType::prism, 3);
    GeometryType pyramid ( GeometryType::pyramid, 3);

    // 1D Quadratures
    Quadrature<double, 1> quadCube1d1(cube1, 1);
    Quadrature<double, 1> quadLine1d1(line, 1);
    Quadrature<double, 1> quadSimplex1d1(simplex1, 1);
    Quadrature<double, 1> quadLine1d3(line, 3);
    Quadrature<double, 1> quadLine1d10(line, 10);

    // 2D Quadratures
    // - Cubes
    Quadrature<double, 2> quadCube2d1(cube2, 1);
    Quadrature<double, 2> quadQuad2d1(quadrilateral, 1);
    Quadrature<double, 2> quadCube2d3(cube2, 3);
    Quadrature<double, 2> quadCube2d9(cube2, 9);
    // - Simplices
    Quadrature<double, 2> quadSimplex2d1(simplex2, 1);
    Quadrature<double, 2> quadTriangle2d1(triangle, 1);
    Quadrature<double, 2> quadSimplex2d4(simplex2, 4);
    Quadrature<double, 2> quadSimplex2d11(simplex2, 11);

    // 3D Quadratures
    // - Cubes
    Quadrature<double, 3> quadCube3d1(cube3, 1);
    Quadrature<double, 3> quadHexa3d1(hexahedron, 1);
    Quadrature<double, 3> quadCube3d2(cube3, 2);
    Quadrature<double, 3> quadCube3d11(cube3, 11);
    // - Simplices
    Quadrature<double, 3> quadSimplex3d1(simplex3, 1);
    Quadrature<double, 3> quadTetra3d1(tetrahedron, 1);
    Quadrature<double, 3> quadSimplex3d4(simplex3, 4);
    Quadrature<double, 3> quadSimplex3d7(simplex3, 7);
    // - Prism and pyramids
    Quadrature<double, 3> quadPrism3d1(prism, 1);
    Quadrature<double, 3> quadPyramid3d1(pyramid, 1);

    /*
    // FixedOrder quads
    typedef FieldVector<double, 2> Coord2;
    typedef FieldVector<double, 3> Coord3;
    FixedOrderQuad<double, Coord2, 1> fixedSimplex2d1(simplex);
    FixedOrderQuad<double, Coord2, 4> fixedSimplex2d4(simplex);
    FixedOrderQuad<double, Coord2, 11> fixedSimplex2d11(simplex);

    FixedOrderQuad<double, Coord3, 1> fixedSimplex3d1(simplex);
    FixedOrderQuad<double, Coord3, 4> fixedSimplex3d4(simplex);
    FixedOrderQuad<double, Coord3, 7> fixedSimplex3d7(simplex);

    //- Test fixed order for simplex
    #ifndef HAVE_ALBERTA
    // if Alberta is present, the new quadratures use the Alberta quadratures
    // while the old ones stick to the UG quadratures
    fixedOrderComparisonExec(quadSimplex2d1, fixedSimplex2d1);
    fixedOrderComparisonExec(quadSimplex2d4, fixedSimplex2d4);
    fixedOrderComparisonExec(quadSimplex2d11, fixedSimplex2d11);
    fixedOrderComparisonExec(quadSimplex3d1, fixedSimplex3d1);
    fixedOrderComparisonExec(quadSimplex3d4, fixedSimplex3d4);
    fixedOrderComparisonExec(quadSimplex3d7, fixedSimplex3d7);
    #endif
    */

    //- Test weight summation for all
    weightSummationExec(quadCube1d1);
    weightSummationExec(quadLine1d1);
    weightSummationExec(quadSimplex1d1);
    weightSummationExec(quadLine1d3);
    weightSummationExec(quadLine1d10);
    weightSummationExec(quadCube2d1);
    weightSummationExec(quadQuad2d1);
    weightSummationExec(quadCube2d3);
    weightSummationExec(quadCube2d9);
    weightSummationExec(quadSimplex2d1);
    weightSummationExec(quadTriangle2d1);
    weightSummationExec(quadSimplex2d4);
    weightSummationExec(quadSimplex2d11);
    weightSummationExec(quadCube3d1);
    weightSummationExec(quadTetra3d1);
    weightSummationExec(quadSimplex3d4);
    weightSummationExec(quadSimplex3d7);
    weightSummationExec(quadPrism3d1);
    weightSummationExec(quadPyramid3d1);

    //- Simple integration test for simplex and cube
    integrationExec(quadCube1d1);
    integrationExec(quadLine1d10);
    integrationExec(quadCube2d3);
    integrationExec(quadCube2d9);
    integrationExec(quadSimplex2d4);
    integrationExec(quadSimplex2d11);
    integrationExec(quadCube3d11);
    integrationExec(quadSimplex3d4);
    integrationExec(quadSimplex3d7);

    //- Test same geometry type
    sameGeometryExec(quadCube1d1, quadLine1d1);
    sameGeometryExec(quadCube2d1, quadQuad2d1);
    sameGeometryExec(quadCube3d1, quadHexa3d1);
    sameGeometryExec(quadSimplex1d1, quadLine1d1);
    sameGeometryExec(quadSimplex2d1, quadTriangle2d1);
    sameGeometryExec(quadSimplex3d1, quadTetra3d1);

    //- Test order
    orderExec(quadCube1d1, 1);
    orderExec(quadLine1d1, 1);
    orderExec(quadSimplex1d1, 1);
    orderExec(quadLine1d3, 3);
    orderExec(quadLine1d10, 10);
    orderExec(quadCube2d1, 1);
    orderExec(quadQuad2d1, 1);
    orderExec(quadCube2d3, 3);
    orderExec(quadCube2d9, 9);
    orderExec(quadSimplex2d1, 1);
    orderExec(quadTriangle2d1, 1);
    orderExec(quadSimplex2d4, 4);
    orderExec(quadSimplex2d11, 11);
    orderExec(quadCube3d1, 1);
    orderExec(quadTetra3d1, 1);
    orderExec(quadSimplex3d4, 4);
    orderExec(quadSimplex3d7, 7); // This one fails, but it's because of ugquadratures.cc

    //- Test uniqueness of indices
    indicesTest();

    std::set<int> ids;
    doTest(ids.insert(quadCube1d1.id()).second);
    doTest(ids.insert(quadLine1d3.id()).second);
    doTest(ids.insert(quadLine1d10.id()).second);
    doTest(ids.insert(quadCube2d1.id()).second);
    doTest(ids.insert(quadCube2d3.id()).second);
    doTest(ids.insert(quadCube2d9.id()).second);
    doTest(ids.insert(quadSimplex2d1.id()).second);
    doTest(ids.insert(quadSimplex2d4.id()).second);
    doTest(ids.insert(quadSimplex2d11.id()).second);
    doTest(ids.insert(quadCube3d1.id()).second);
    doTest(ids.insert(quadCube3d2.id()).second);
    doTest(ids.insert(quadCube3d11.id()).second);
    doTest(ids.insert(quadSimplex3d1.id()).second);
    doTest(ids.insert(quadSimplex3d4.id()).second);
    doTest(ids.insert(quadSimplex3d7.id()).second);
    doTest(ids.insert(quadPrism3d1.id()).second);
    doTest(ids.insert(quadPyramid3d1.id()).second);

  }

  template <class Quad, class Fixed>
  void Quad_Test::fixedOrderComparisonExec(Quad& quad, Fixed& fixed)
  {
    for (int i = 0; i < quad.nop(); ++i) {
      for (int j = 0; j < Quad::dimension; ++j) {
        doTest(quad.point(i)[j], fixed.point(i)[j]);
      }
      doTest(quad.weight(i), fixed.weight(i));
    }
  }

  template <class Quad>
  void Quad_Test::weightSummationExec(Quad& quad)
  {
    double sum = 0.;
    for (int i = 0; i < quad.nop(); ++i) {
      sum += quad.weight(i);
    }

    const GeometryType geom = quad.geometry();
    if ( geom.isCube())
    {
      doTest(sum, 1.0);
    }
    else if ( geom.isSimplex())
    {
      switch(Quad::dimension) {
      case 1:
        doTest(sum, 1.0);
        break;
      case 2:
        doTest(sum, 0.5);
        break;
      case 3:
        doTestTol(sum, 0.16666666666666666667, 1e-04);
        break;
      }
    }
    else if( geom.isPrism() )
    {
      doTest(sum, 0.5);
    }
    else if ( geom.isPyramid() )
    {
      doTestTol(sum, 0.3333333333333333, 1e-04);
    }
    else
    {
      DUNE_THROW(NotImplemented, "What geometry type ey?");
    }
  }

  template <class Quad>
  void Quad_Test::integrationExec(Quad& quad)
  {
    Func f;
    double result = 0.;
    for (int i = 0; i < quad.nop(); ++i) {
      result += f(quad.point(i))*quad.weight(i);
    }

    switch (Quad::dimension) {
    case 1:
      doTest(result, -0.25);
      break;
    case 2:
      if (quad.geometry().isCube() ) {
        doTest(result, -0.2708333333);
      }
      else {
        doTest(result, -0.0333333333);
      }
      break;
    case 3:
      if (quad.geometry().isCube()) {
        doTest(result, -0.03854166666);
      }
      else {
        doTestTol(result, -0.00551388888, 1e-05);
      }
      break;
    default:
      _fail("should not get here");
    }
  }

  template <class Quad>
  void Quad_Test::orderExec(Quad& quad, int order)
  {
    doTest(quad.order() >= order);
  }

  template <class Quad>
  void Quad_Test::sameGeometryExec(Quad& quad1, Quad& quad2)
  {
    doTest(quad1.id() == quad2.id());
  }

  void Quad_Test::indicesTest()
  {
    const GeometryType line ( GeometryType::cube, 1 );
    size_t id;
    {
      Quadrature<double, 1> quadTemp(line, 5);
      id = quadTemp.id();
    }
    Quadrature<double, 1> quadTemp(line, 5);
    doTest(quadTemp.id() == id);
  }
  } // end namespace Fem
} // end namespace Dune
