#include <config.h>

#define USE_TWISTFREE_MAPPER

#define TEST_SECOND_ORDER
#define TEST_THIRD_ORDER

#ifdef ALUGRID_SIMPLEX
#if GRIDDIM > 2
#undef TEST_THIRD_ORDER
#endif
#endif

#include <dune/fem/misc/suite.hh>
using namespace Dune;
using namespace Fem;

#include "basisfunctiontest.cc"
#include "mappertest.cc"

int main( int argc, char **argv )
{
  MPIManager::initialize( argc, argv );

  std::stringstream gridFile;
  gridFile << GRIDDIM<<"dgrid.dgf";

  Suite suite("Basisfunction tests");
  suite.addTest( new LagrangeBasis_Test(gridFile.str()) );
  suite.addTest( new LagrangeMapper_Test< GridSelector::GridType >( gridFile.str() ) );

  suite.run();
  suite.report();
}
