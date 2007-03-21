#include <dune/fem/misc/suite.hh>

#include <dune/config.h>

using namespace Dune;

#include "basefunctiontest.cc"

int main() 
{
  std :: string gridFile( "../../test/2dgrid.dgf" );
//  std :: string gridFile( "../../../examples/poisson/3dgrid.al" );
  
  Suite suite("Basefunction tests");
  suite.addTest(new LagrangeBase_Test(gridFile));

  suite.run();
  suite.report();
}
