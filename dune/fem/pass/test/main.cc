#include <iostream>
#include <string>

#include <config.h>

#ifndef GRIDDIM 
#define GRIDDIM 3
#endif
const int dimension = GRIDDIM;

#include "suite.hh"
#include <dune/fem/misc/test.hh>

#include "pass_test.hh"
#include "helper_test.hh"

using namespace Dune;

int main() {
  Suite passSuite("Test suite for pass implementation");

  std::stringstream filename;
  filename << "dgellipt/grid" << dimension << "d.dgf";
  passSuite.addTest(new Pass_Test(filename.str()));
  passSuite.addTest(new PassHelper_Test(filename.str()));
  
  passSuite.run();
  passSuite.report();
}

