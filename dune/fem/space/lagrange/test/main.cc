#include <config.h>

#define USE_TWISTFREE_MAPPER

#define TEST_SECOND_ORDER
#define TEST_THIRD_ORDER

#include <dune/fem/misc/suite.hh>
using namespace Dune;
using namespace Fem;

#include "basisfunctiontest.cc"
#include "mappertest.cc"

int main( int argc, char **argv ) 
{
  MPIManager::initialize( argc, argv );
  //std :: string gridFile( "../../../examples/elliptic/square.dgf" );
  //std :: string gridFile( "onesimplex.dgf" );
  //std :: string gridFile( "cube.dgf" );
  //std :: string gridFile( "../../test/2dgrid.dgf" );
  std :: string gridFile( "2dgrid.dgf" );
  //std :: string gridFile( "../../../examples/poisson/3dgrid.al" );
  
  Suite suite("Basisfunction tests");
  suite.addTest( new LagrangeBasis_Test(gridFile));
  suite.addTest( new LagrangeMapper_Test< GridSelector::GridType >( gridFile ) );

  suite.run();
  suite.report();
}
