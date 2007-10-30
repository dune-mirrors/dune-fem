#include <iostream>

#include <config.h>

#include "../../misc/test.hh"
#include "../../misc/suite.hh"

// if grid dim is not defined simply use ALBERTA_DIM 
#ifndef GRIDDIM 
#define ALBERTA_DIM GRIDDIM
#endif

#include "cachequad_test.hh"
#include "quad_test.hh"

int main() {

  std::stringstream albertaGridFile;
  albertaGridFile << "../../../macrogrids/AlbertaGrid/" << GRIDDIM << "dgrid.al";
  std::string aluGridHexaFile("../../../macrogrids/ALU3dGrid/cube.hexa");
  std::string aluGridTetraFile("../../../macrogrids/ALU3dGrid/macro.small");

  Dune::Suite quadSuite("Tests for new quadratures");
  
  quadSuite.addTest(new Dune::Quad_Test());
  quadSuite.addTest(new Dune::CachingQuadrature_Test(albertaGridFile.str(),
                                                     aluGridHexaFile,
                                    						     aluGridTetraFile));

  quadSuite.run();
  quadSuite.report();
}
