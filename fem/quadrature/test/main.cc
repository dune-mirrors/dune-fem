#include <iostream>

#include <config.h>

#include "include_cc.hh"

#include "test.hh"
#include "suite.hh"
#include "quad_test.hh"

using namespace Dune;

int main() {
  Suite quadSuite("Tests for new quadratures");
  
  quadSuite.addTest(new Quad_Test());

  quadSuite.run();
  quadSuite.report();
}
