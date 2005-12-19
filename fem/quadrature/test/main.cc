#include <iostream>

//#include <config.h>

#include "../../misc/test.hh"
#include "../../misc/suite.hh"

#include "cachequad_test.hh"
#include "quad_test.hh"

#include "include_cc.hh"

int main() {
  std::string albertaGridFile("../../macrogrids/AlbertaGrid/2dmacro.al");
  std::string aluGridHexaFile("../../macrogrids/ALU3dGrid/cube.hexa");
  std::string aluGridTetraFile("../../macrogrids/ALU3dGrid/macro.small");

  Dune::Suite quadSuite("Tests for new quadratures");
  
  quadSuite.addTest(new Dune::Quad_Test());
  quadSuite.addTest(new Dune::CachingQuadrature_Test(albertaGridFile,
                                                     aluGridHexaFile,
						     aluGridTetraFile));

  quadSuite.run();
  quadSuite.report();
}
