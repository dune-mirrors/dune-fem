#include <config.h>

#include "twist_test.hh"

#include "../../quadrature.hh"

namespace Dune {
  
  void TwistProvider_Test::run() 
  {
    lineTest();
    triangleTest();
    quadrilateralTest();
    nonSymmetricTest();
  }

  void TwistProvider_Test::lineTest() 
  {

    GeometryType line(GeometryType::simplex,1);


    typedef FieldVector<double, 1> CoordinateType;
    typedef TwistProvider<double, 1>::TwistStorageType TwistStorageType;
    typedef TwistStorageType::PointVectorType PointVectorType;
    typedef TwistStorageType::MapperType MapperType;

    CoordinateType begin(0.1);
    CoordinateType end(0.9);

    TestQuadrature<double, 1> quadImp(line, 0);
    quadImp.newQuadraturePoint(begin, 0.5);
    quadImp.newQuadraturePoint(end, 0.5);
    
    Quadrature<double, 1> quad(quadImp);
 
    const TwistStorageType& storage =
      TwistProvider<double, 1>::getTwistStorage(quad.ipList());

    // Test minTwist, maxTwist
    _test(storage.minTwist() == 0);
    _test(storage.maxTwist() == 2);
    
    // Test points
    const PointVectorType& points = storage.getPoints();
    _floatTest(points[0][0], begin[0]);
    _floatTest(points[1][0], end[0]);

    // Test mapping
    const MapperType& m0 = storage.getMapper(0);
    _test(m0[0] == 0);
    _test(m0[1] == 1);
    
    const MapperType& m1 = storage.getMapper(1);
    _test(m1[0] == 1);
    _test(m1[1] == 0);

  }

  void TwistProvider_Test::triangleTest() 
  {
    GeometryType simplex(GeometryType::simplex,2);
    typedef FieldVector<double, 2> CoordinateType;
    typedef TwistProvider<double, 2>::TwistStorageType TwistStorageType;
    typedef TwistStorageType::PointVectorType PointVectorType;
    typedef TwistStorageType::MapperType MapperType;

    CoordinateType p0(0.0);
    CoordinateType p1(0.0); p1[0] = 1.0;
    CoordinateType p2(0.0); p2[1] = 1.0;

    const double oneThird = 1.0/3.0;

    TestQuadrature<double, 2> quadImp(simplex, 0);
    quadImp.newQuadraturePoint(p0, oneThird);
    quadImp.newQuadraturePoint(p1, oneThird);
    quadImp.newQuadraturePoint(p2, oneThird);

    Quadrature<double, 2> quad(quadImp);
  
    const TwistStorageType& storage =
      TwistProvider<double, 2>::getTwistStorage(quad.ipList());

    // Test minTwist, maxTwist
    _test(storage.minTwist() == -3);
    _test(storage.maxTwist() == 3);

    // Test points
    const PointVectorType& points = storage.getPoints();
    _floatTest(points[0][0], p0[0]);
    _floatTest(points[0][1], p0[1]);
    _floatTest(points[1][0], p1[0]);
    _floatTest(points[1][1], p1[1]);
    _floatTest(points[2][0], p2[0]);
    _floatTest(points[2][1], p2[1]);

    // Test twists
    const MapperType& m_3 = storage.getMapper(-3);
    _test(m_3[0] == 2);
    _test(m_3[1] == 1);
    _test(m_3[2] == 0);
    //std::cout << "(" << m_3[0] << ", " << m_3[1] << ", " << m_3[2] << ")\n";

    const MapperType& m_2 = storage.getMapper(-2);
    _test(m_2[0] == 1);
    _test(m_2[1] == 0);
    _test(m_2[2] == 2);
    //std::cout << "(" << m_2[0] << ", " << m_2[1] << ", " << m_2[2] << ")\n";

    const MapperType& m_1 = storage.getMapper(-1);
    _test(m_1[0] == 0);
    _test(m_1[1] == 2);
    _test(m_1[2] == 1);
    //std::cout << "(" << m_1[0] << ", " << m_1[1] << ", " << m_1[2] << ")\n";

    const MapperType& m0 = storage.getMapper(0);
    _test(m0[0] == 0);
    _test(m0[1] == 1);
    _test(m0[2] == 2);
    //std::cout << "(" << m0[0] << ", " << m0[1] << ", " << m0[2] << ")\n";

    const MapperType& m1 = storage.getMapper(1);
    _test(m1[0] == 1);
    _test(m1[1] == 2);
    _test(m1[2] == 0);
    //std::cout << "(" << m1[0] << ", " << m1[1] << ", " << m1[2] << ")\n";

    const MapperType& m2 = storage.getMapper(2);
    _test(m2[0] == 2);
    _test(m2[1] == 0);
    _test(m2[2] == 1);
    //std::cout << "(" << m2[0] << ", " << m2[1] << ", " << m2[2] << ")\n";
  }

  void TwistProvider_Test::quadrilateralTest() 
  {
    GeometryType cube(GeometryType::cube,2);
    typedef FieldVector<double, 2> CoordinateType;
    typedef TwistProvider<double, 2>::TwistStorageType TwistStorageType;
    typedef TwistStorageType::PointVectorType PointVectorType;
    typedef TwistStorageType::MapperType MapperType;

    CoordinateType p0(0.0);
    CoordinateType p1(0.0); p1[0] = 1.0;
    CoordinateType p2(0.0); p2[1] = 1.0;
    CoordinateType p3(1.0); 

    TestQuadrature<double, 2> quadImp(cube, 0);
    quadImp.newQuadraturePoint(p0, 0.25);
    quadImp.newQuadraturePoint(p1, 0.25);
    quadImp.newQuadraturePoint(p2, 0.25);
    quadImp.newQuadraturePoint(p3, 0.25);

    Quadrature<double, 2> quad(quadImp);
  
    const TwistStorageType& storage =
      TwistProvider<double, 2>::getTwistStorage(quad.ipList());

    // Test minTwist, maxTwist
    _test(storage.minTwist() == -4);
    _test(storage.maxTwist() == 4);

    // Test points
    const PointVectorType& points = storage.getPoints();
    _floatTest(points[0][0], p0[0]);
    _floatTest(points[0][1], p0[1]);
    _floatTest(points[1][0], p1[0]);
    _floatTest(points[1][1], p1[1]);
    _floatTest(points[2][0], p2[0]);
    _floatTest(points[2][1], p2[1]);
    _floatTest(points[3][0], p3[0]);
    _floatTest(points[3][1], p3[1]);
    
    // Test twists
    const MapperType& m_4 = storage.getMapper(-4);
    _test(m_4[0] == 2);
    _test(m_4[1] == 3);
    _test(m_4[2] == 0);
    _test(m_4[3] == 1);
    //std::cout << "(" << m_4[0] << ", " << m_4[1] << ", " << m_4[2] << ", " << m_4[3] << ")\n";

    const MapperType& m_3 = storage.getMapper(-3);
    _test(m_3[0] == 3);
    _test(m_3[1] == 1);
    _test(m_3[2] == 2);
    _test(m_3[3] == 0);
    //std::cout << "(" << m_3[0] << ", " << m_3[1] << ", " << m_3[2] << ", " << m_3[3] << ")\n";

    const MapperType& m_2 = storage.getMapper(-2);
    _test(m_2[0] == 1);
    _test(m_2[1] == 0);
    _test(m_2[2] == 3);
    _test(m_2[3] == 2);
    //std::cout << "(" << m_2[0] << ", " << m_2[1] << ", " << m_2[2] << ", " << m_2[3] << ")\n";

    const MapperType& m_1 = storage.getMapper(-1);
    _test(m_1[0] == 0);
    _test(m_1[1] == 2);
    _test(m_1[2] == 1);
    _test(m_1[3] == 3);
    // std::cout << "(" << m_1[0] << ", " << m_1[1] << ", " << m_1[2] << ", " << m_1[3] << ")\n";

    const MapperType& m0 = storage.getMapper(0);
    _test(m0[0] == 0);
    _test(m0[1] == 1);
    _test(m0[2] == 2);
    _test(m0[3] == 3);
    //std::cout << "(" << m0[0] << ", " << m0[1] << ", " << m0[2] << ", " << m0[3] << ")\n";

    const MapperType& m1 = storage.getMapper(1);
    _test(m1[0] == 1);
    _test(m1[1] == 3);
    _test(m1[2] == 0);
    _test(m1[3] == 2);
    // std::cout << "(" << m1[0] << ", " << m1[1] << ", " << m1[2] << ", " << m1[3] << ")\n";

    const MapperType& m2 = storage.getMapper(2);
    _test(m2[0] == 3);
    _test(m2[1] == 2);
    _test(m2[2] == 1);
    _test(m2[3] == 0);
    //std::cout << "(" << m2[0] << ", " << m2[1] << ", " << m2[2] << ", " << m2[3] << ")\n";

    const MapperType& m3 = storage.getMapper(3);
    _test(m3[0] == 2);
    _test(m3[1] == 0);
    _test(m3[2] == 3);
    _test(m3[3] == 1);
    //std::cout << "(" << m3[0] << ", " << m3[1] << ", " << m3[2] << ", " << m3[3] << ")\n";
  }

  void TwistProvider_Test::nonSymmetricTest() 
  {
    typedef FieldVector<double, 1> CoordinateType;
    typedef TwistProvider<double, 1>::TwistStorageType TwistStorageType;
    typedef TwistStorageType::PointVectorType PointVectorType;
    typedef TwistStorageType::MapperType MapperType;

    CoordinateType begin(0.1);
    CoordinateType end(0.9);
    
    GeometryType line(GeometryType::simplex,1);
    TestQuadrature<double, 1> quadImp(line, 0);
    quadImp.newQuadraturePoint(begin, 0.5);
    quadImp.newQuadraturePoint(end, 0.5);
    
    Quadrature<double, 1> quad(quadImp);
    
    const TwistStorageType& storage =
      TwistProvider<double, 1>::getTwistStorage(quad.ipList());
    
    // Test points
    const PointVectorType& points = storage.getPoints();
    _test(points.size() == 2);
    _floatTest(points[0][0], begin[0]);
    _floatTest(points[1][0], end[0]);
    
    // Test mapping
    const MapperType& m0 = storage.getMapper(0);
    _test(m0[0] == 0);
    
    const MapperType& m1 = storage.getMapper(1);
    _test(m1[0] == 1);
  
  }
}
