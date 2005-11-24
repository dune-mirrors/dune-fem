#include "quad_test.hh"

namespace Dune {

  void Quad_Test::run() {
    _fail("No tests so far");

    Quadrature<double, 1> quad1(cube, 1);
    Quadrature<double, 2> quad2(cube, 1);
    Quadrature<double, 3> quad3(cube, 1);
  }

}
