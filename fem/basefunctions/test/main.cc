#include "../../misc/suite.hh"

#include <dune/config.h>

#include "include.cc"
#include "storage_test.hh"

using namespace Dune;

int main() 
{
  std::string gridFile("../../macrogrids/ALU3dGrid/cube.hexa");
  //std::string gridFile("../../macrogrids/AlbertaGrid/2dmacro.al");
  Suite suite("Basefunction tests");

  suite.addTest(new Storage_Test(gridFile));

  suite.run();
  suite.report();
}
