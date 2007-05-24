#include <dune/fem/misc/suite.hh>

#include <dune/config.h>

using namespace Dune;

#include "basefunctiontest.cc"
#include "mappertest.cc"

int main() 
{
  //std :: string gridFile( "../../../examples/elliptic/square.dgf" );
  //std :: string gridFile( "onesimplex.dgf" );
  //std :: string gridFile( "cube.dgf" );
  std :: string gridFile( "../../test/2dgrid.dgf" );
  //std :: string gridFile( "../../../examples/poisson/3dgrid.al" );
  
  Suite suite("Basefunction tests");
  suite.addTest(new LagrangeBase_Test(gridFile));
  suite.addTest( new LagrangeMapper_Test( gridFile ) );

  suite.run();
  suite.report();
}
