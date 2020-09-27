#include <config.h>

#include "point_test.hh"

namespace Dune {
  namespace Fem {

  void PointProvider_Test::run()
  {
    codim0Test();
    sameOutputTest();
    // TODO: fix transformationTest
    // transformationTest();
  }


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



  void PointProvider_Test::codim0Test()
  {
    const int dim = 2;
    typedef PointProvider<double,dim,0>::GlobalPointVectorType PointVectorType;

    GeometryType simplex = GeometryTypes::simplex( dim );

    Quadrature<double, dim> quad(simplex, dim);

    PointProvider<double, dim, 0>::registerQuadrature(quad.ipList());

    const PointVectorType& points =
      PointProvider<double, dim, 0>::getPoints(quad.id(), simplex);

    doTest( (int) points.size() == quad.nop() );
    for (size_t i = 0; i < points.size(); ++i)
    {
      for (int j = 0; j < dim; ++j) {
        doTest(points[i][j], quad.point(i)[j]);
      }
    }
  }

  void PointProvider_Test::sameOutputTest()
  {
    const int dim = 2;
    const int codim = 1;
    typedef PointProvider<double, dim, codim> PointProviderType;
    typedef PointProviderType::GlobalPointVectorType PointVectorType;
    typedef FieldVector<double, dim-codim> LocalPointType;

    LocalPointType first(0.5);
    LocalPointType second(0.6);

    GeometryType simplex = GeometryTypes::simplex( dim );
    GeometryType line = GeometryTypes::simplex( dim-1 );

    TestQuadrature<double, 1> quadImp(line, 0);
    quadImp.newQuadraturePoint(first, 0.5);
    quadImp.newQuadraturePoint(second, 0.5);

    Quadrature<double, 1> quad1(quadImp);
    Quadrature<double, 1> quad2(quadImp);

    PointProviderType::getMappers(quad1.ipList(), simplex);
    PointProviderType::getMappers(quad2.ipList(), simplex);

    const PointVectorType& p1 =
      PointProviderType::getPoints(quad1.id(), simplex);
    const PointVectorType& p2 =
      PointProviderType::getPoints(quad2.id(), simplex);

    doTest(p1.size() == p2.size());
    for (size_t i = 0; i < p1.size(); ++i) {
      for (int j = 0; j < dim; ++j) {
        doTest(p1[i][j], p2[i][j]);
      }
    }
  }

  void PointProvider_Test::transformationTest()
  {
    GeometryType quadrilateral = GeometryTypes::quadrilateral;
    GeometryType triangle = GeometryTypes::triangle;
    GeometryType lineC = GeometryTypes::cube( 1 );

    GeometryType lineS = GeometryTypes::simplex( 1 );

    typedef PointProvider<double, 2, 1> PointProvider1Type;
    typedef PointProvider<double, 3, 1> PointProvider2Type;
    typedef FieldVector<double, 1> Point1Type;
    typedef FieldVector<double, 2> Point2Type;
    typedef FieldVector<double, 3> Point3Type;

    Point2Type tmp2(0.5);
    Point3Type tmp3(0.5);

    const double oneThird = 1.0/3.0;

    // Hexa test
    Point2Type ph(0.5);

    TestQuadrature<double, 2> quadImpQuad(quadrilateral, 0);
    quadImpQuad.newQuadraturePoint(ph, 1.0);
    Quadrature<double, 2> quadQuad(quadImpQuad);

    const auto& mvh =
      PointProvider2Type::getMappers(quadQuad.ipList() , quadrilateral);
    for (size_t i = 0; i < mvh.first.size(); ++i) {
      std::cout << mvh.first[i][0] << ", ";
    }
    std::cout << std::endl;


    const PointProvider2Type::GlobalPointVectorType& ptsHexa =
      PointProvider2Type::getPoints(quadQuad.id(), quadrilateral);

    doTest(ptsHexa.size() == 6);
    tmp3[0] = 0.;
    doTest(findPoint(tmp3, ptsHexa));
    tmp3[0] = 0.5; tmp3[1] = 0.;
    doTest(findPoint(tmp3, ptsHexa));
    tmp3[1] = 0.5; tmp3[2] = 0.;
    doTest(findPoint(tmp3, ptsHexa));
    tmp3[2] = 1.0;
    doTest(findPoint(tmp3, ptsHexa));
    tmp3[2] = 0.5; tmp3[1] = 1.0;
    doTest(findPoint(tmp3, ptsHexa));
    tmp3[1] = 0.5; tmp3[0] = 1.0;
    doTest(findPoint(tmp3, ptsHexa));

    // Tetra test
    Point2Type pt(oneThird);
    TestQuadrature<double, 2> quadImpTri(triangle, 0);
    quadImpTri.newQuadraturePoint(pt, 1.0);
    Quadrature<double, 2> quadTri(quadImpTri);

    PointProvider2Type::getMappers(quadTri.ipList() , triangle);
    const PointProvider2Type::GlobalPointVectorType& ptsTetra =
      PointProvider2Type::getPoints(quadTri.id(), triangle);

    doTest(ptsTetra.size() == 4);
    tmp3 = oneThird;
    doTest(findPoint(tmp3, ptsTetra));
    tmp3[0] = 0.;
    doTest(findPoint(tmp3, ptsTetra));
    tmp3[1] = 0.; tmp3[0] = oneThird;
    doTest(findPoint(tmp3, ptsTetra));
    tmp3[2] = 0.; tmp3[1] = oneThird;
    doTest(findPoint(tmp3, ptsTetra));

     // Quadrilateral test
    Point1Type pl(0.5);
    TestQuadrature<double, 1> quadImpLine(lineC, 0);
    quadImpLine.newQuadraturePoint(pl, 1.0);
    Quadrature<double, 1> quadLine(quadImpLine);

    PointProvider1Type::getMappers(quadLine.ipList(), quadrilateral );
    const PointProvider1Type::GlobalPointVectorType& ptsQuad =
      PointProvider1Type::getPoints(quadLine.id(), quadrilateral);

    doTest(ptsQuad.size() == 4);
    tmp2[0] = 0.;
    doTest(findPoint(tmp2, ptsQuad));
    tmp2[1] = 0.; tmp2[0] = 0.5;
    doTest(findPoint(tmp2, ptsQuad));
    tmp2[0] = 1.; tmp2[1] = 0.5;
    doTest(findPoint(tmp2, ptsQuad));
    tmp2[1] = 1.; tmp2[0] = 0.5;
    doTest(findPoint(tmp2, ptsQuad));

    // Triangle test
    TestQuadrature<double, 1> quadImpLine2(lineS, 0);
    quadImpLine2.newQuadraturePoint(pl, 1.0);

    Quadrature<double, 1> quadLine2(quadImpLine2);
    PointProvider1Type::getMappers(quadLine2.ipList() , triangle);

    const PointProvider1Type::GlobalPointVectorType& ptsTri =
      PointProvider1Type::getPoints(quadLine2.id(), triangle);

    doTest(ptsTri.size() == 3);
    tmp2 = 0.5;
    doTest(findPoint(tmp2, ptsTri));
    tmp2[0] = 0.;
    doTest(findPoint(tmp2, ptsTri));
    tmp2[0] = 0.5; tmp2[1] = 0;
    doTest(findPoint(tmp2, ptsTri));
  }

  template <class PointType>
  bool PointProvider_Test::findPoint(const PointType& p,
                                      const std::vector<PointType>& vec)
  {
    const double eps = 1.0e-5;

    for (size_t i = 0; i < vec.size(); ++i) {
      bool same = false;
      for (size_t j = 0; j < p.dim(); ++j) {
        same |= (p[j] > vec[i][j]-eps && p[j] < vec[i][j]+eps);
      }

      if (same) {
        return true;
      }
    }
    return false;
  }

  } // end namespace Fem
} // end namespace Dune
