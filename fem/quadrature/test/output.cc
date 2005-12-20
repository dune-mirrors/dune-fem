#include <iostream>
#include <iomanip>

#include <dune/quadrature/fixedorder.hh>
#include <dune/common/misc.hh>

const int maxOrder = 20;
const int precision = 50;

using namespace Dune;

template <class Quad>
void print(Quad& quad) {
  std::cout.precision(precision);
  std::cout << "Order:\n" << quad.order() << std::endl;
  std::cout << "Number of points:\n" << quad.nop() << std::endl;
  std::cout << "Points:\n{\n";
  for (int i = 0; i < quad.nop(); ++i) {
    std::cout << quad.point(i) << ",\n";
  }
  std::cout << "}" << std::endl;
  std::cout << "Weights:\n{\n";
  for (int i = 0; i < quad.nop(); ++i) {
    std::cout << quad.weight(i) << "\n";
  }
  std::cout << "}" << std::endl;
   
}

template <int dim>
void outputQuadrature(GeometryType geo) {
  for (int i = 1; i < maxOrder; ++i) {
    QuadratureOld<double, FieldVector<double, 2> > quad(0, geo, i);
    if (i > quad.order()) break;
    print(quad);
    std::cout << std::endl;
  }
      
}

int main() {
    
  GeometryType geo;

  geo = triangle;
  std::cout << "Triangle output\n"
            << "---------------\n";
  try {
    outputQuadrature<2>(geo);
  } catch(...) {}

  geo = tetrahedron;
  std::cout << "Tetrahedron output\n"
            << "---------------\n";
  try {
    outputQuadrature<3>(geo);
  } catch(...) {}

  geo = prism;
  std::cout << "Prism output\n"
            << "---------------\n";
  try {
    outputQuadrature<3>(geo);
  } catch(...) {}

  geo = pyramid;
  std::cout << "Pyramid output\n"
            << "---------------\n";
  try {
    outputQuadrature<3>(geo);
  } catch(...) {}
}
