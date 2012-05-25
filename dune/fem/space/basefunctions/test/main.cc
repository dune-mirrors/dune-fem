#include <config.h>

#include <dune/fem/misc/suite.hh>
#include "storage_test.hh"

using namespace Dune;

int main() 
{
  //std::string gridFile("../../macrogrids/ALU3dGrid/cube.hexa");
  std::string gridFile("2dgrid.dgf");
  Suite suite("Basefunction tests");

  suite.addTest(new Storage_Test(gridFile));

  suite.run();
  suite.report();
}
