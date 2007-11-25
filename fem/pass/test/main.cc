#include <iostream>

#include <dune/config.h>

#include "suite.hh"
#include <dune/fem/misc/test.hh>
#include "pass_test.hh"
#include "helper_test.hh"

using namespace Dune;

int main() {
  Suite passSuite("Test suite for pass implementation");

  passSuite.addTest(new Pass_Test("macro.small"));
  passSuite.addTest(new PassHelper_Test("macro.small"));
  
  passSuite.run();
  passSuite.report();
}

