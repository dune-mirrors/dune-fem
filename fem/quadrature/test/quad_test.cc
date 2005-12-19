#include <config.h>

#include "quad_test.hh"
#include <set>

namespace Dune {

  void Quad_Test::run() {
    // 1D Quadratures
    Quadrature<double, 1> quadCube1d1(cube, 1);
    Quadrature<double, 1> quadLine1d1(line, 1);
    Quadrature<double, 1> quadSimplex1d1(simplex, 1);
    Quadrature<double, 1> quadLine1d3(line, 3);
    Quadrature<double, 1> quadLine1d10(line, 10);

    // 2D Quadratures
    // - Cubes
    Quadrature<double, 2> quadCube2d1(cube, 1);
    Quadrature<double, 2> quadQuad2d1(quadrilateral, 1);
    Quadrature<double, 2> quadCube2d3(cube, 3);
    Quadrature<double, 2> quadCube2d9(cube, 9);
    // - Simplices
    Quadrature<double, 2> quadSimplex2d1(simplex, 1);
    Quadrature<double, 2> quadTriangle2d1(triangle, 1);
    Quadrature<double, 2> quadSimplex2d4(simplex, 4);
    Quadrature<double, 2> quadSimplex2d11(simplex, 11);

    // 3D Quadratures
    // - Cubes
    Quadrature<double, 3> quadCube3d1(cube, 1);
    Quadrature<double, 3> quadHexa3d1(hexahedron, 1);
    Quadrature<double, 3> quadCube3d2(cube, 2);
    Quadrature<double, 3> quadCube3d11(cube, 11);
    // - Simplices
    Quadrature<double, 3> quadSimplex3d1(simplex, 1);
    Quadrature<double, 3> quadTetra3d1(tetrahedron, 1);
    Quadrature<double, 3> quadSimplex3d4(simplex, 4);
    Quadrature<double, 3> quadSimplex3d7(simplex, 7);
    // - Prism and pyramids
    Quadrature<double, 3> quadPrism3d1(prism, 1);
    Quadrature<double, 3> quadPyramid3d1(pyramid, 1);

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
    _test(ids.insert(quadCube1d1.id()).second);
    _test(ids.insert(quadLine1d3.id()).second);
    _test(ids.insert(quadLine1d10.id()).second);
    _test(ids.insert(quadCube2d1.id()).second);
    _test(ids.insert(quadCube2d3.id()).second);
    _test(ids.insert(quadCube2d9.id()).second);
    _test(ids.insert(quadSimplex2d1.id()).second);
    _test(ids.insert(quadSimplex2d4.id()).second);
    _test(ids.insert(quadSimplex2d11.id()).second);
    _test(ids.insert(quadCube3d1.id()).second);
    _test(ids.insert(quadCube3d2.id()).second);
    _test(ids.insert(quadCube3d11.id()).second);
    _test(ids.insert(quadSimplex3d1.id()).second);
    _test(ids.insert(quadSimplex3d4.id()).second);
    _test(ids.insert(quadSimplex3d7.id()).second);
    _test(ids.insert(quadPrism3d1.id()).second);
    _test(ids.insert(quadPyramid3d1.id()).second);

  }

  template <class Quad, class Fixed>
  void Quad_Test::fixedOrderComparisonExec(Quad& quad, Fixed& fixed) 
  {
    for (int i = 0; i < quad.nop(); ++i) {
      for (int j = 0; j < Quad::dimension; ++j) {
        _floatTest(quad.point(i)[j], fixed.point(i)[j]);
      }
      _floatTest(quad.weight(i), fixed.weight(i));
    }
  }

  template <class Quad>
  void Quad_Test::weightSummationExec(Quad& quad) 
  {
    double sum = 0.;  
    for (int i = 0; i < quad.nop(); ++i) {
      sum += quad.weight(i);
    }

    switch(quad.geometry()) {
    case line:
      _floatTest(sum, 1.0);
      break;
    case cube:
      _floatTest(sum, 1.0);
      break;
    case simplex:
      switch(Quad::dimension) {
      case 1:
        _floatTest(sum, 1.0);
        break;
      case 2:
        _floatTest(sum, 0.5);
        break;
      case 3:
        _floatTestTol(sum, 0.16666666666666666667, 1e-04);
        break;
      }
      break;
    case triangle:
      _floatTest(sum, 0.5);
      break;
    case quadrilateral:
      _floatTest(sum, 1.0);
      break;
    case tetrahedron:
      _floatTest(sum, 0.16666666666666667);
      break;
    case hexahedron:
      _floatTest(sum, 1.0);
      break;
    case prism:
      _floatTest(sum, 0.5);
      break;
    case pyramid:
      _floatTestTol(sum, 0.3333333333333333, 1e-04);
      break;
    default:
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
      _floatTest(result, -0.25);
      break;
    case 2:
      if (quad.geometry() == cube) {
        _floatTest(result, -0.2708333333);
      } 
      else {
        _floatTest(result, -0.0333333333);
      } 
      break;
    case 3:
      if (quad.geometry() == cube) {
        _floatTest(result, -0.03854166666);
      } 
      else {
        _floatTestTol(result, -0.00551388888, 1e-05);
      } 
      break;
    default:
      _fail("should not get here");
    }
  }

  template <class Quad>
  void Quad_Test::orderExec(Quad& quad, int order) 
  {
    _test(quad.order() >= order);
  }
  
  template <class Quad>
  void Quad_Test::sameGeometryExec(Quad& quad1, Quad& quad2) 
  {
    _test(quad1.id() == quad2.id());
  }

  void Quad_Test::indicesTest() 
  {
    size_t id;
    {
      Quadrature<double, 1> quadTemp(line, 5);
      id = quadTemp.id();
    }
    Quadrature<double, 1> quadTemp(line, 5);
    _test(quadTemp.id() == id);
  }
}
