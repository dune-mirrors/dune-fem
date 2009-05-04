#include <iostream>

#include <config.h>

#include <dune/fem/misc/mpimanager.hh>

#include "../../misc/test.hh"
#include "../../misc/suite.hh"

#include "cachequad_test.hh"
#include "quad_test.hh"

int main ( int argc, char **argv )
{
  Dune::MPIManager::initialize( argc, argv );

  std::stringstream albertaGridFile;
#if HAVE_ALBERTA
  albertaGridFile << "../../../macrogrids/AlbertaGrid/" << ALBERTA_DIM << "dgrid.al";
#endif

  std::string aluGridHexaFile("../../../macrogrids/ALU3dGrid/cube.hexa");
  std::string aluGridTetraFile("../../../macrogrids/ALU3dGrid/macro.small");
  std::string dgf2DGridFile("../../../macrogrids/DGFMacrogrids/examplegrid5.dgf");
  std::string dgf3DGridFile("../../../macrogrids/DGFMacrogrids/examplegrid6.dgf");

  Dune::Suite quadSuite("Tests for new quadratures");
  
  quadSuite.addTest(new Dune::Quad_Test());
  quadSuite.addTest(new Dune::CachingQuadrature_Test(albertaGridFile.str(),
                                                     aluGridHexaFile,
                                    						     aluGridTetraFile,
                                                     dgf2DGridFile,
                                                     dgf3DGridFile));

  quadSuite.run();
  quadSuite.report();
}
