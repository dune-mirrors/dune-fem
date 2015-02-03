#include <config.h>

#define USE_TWISTFREE_MAPPER

#define TEST_SECOND_ORDER
#define TEST_THIRD_ORDER

#ifdef ALUGRID_SIMPLEX
#if GRIDDIM > 2
#undef TEST_THIRD_ORDER
#endif
#endif

#include <dune/fem/misc/mpimanager.hh>

#include "basisfunctiontest.cc"
#include "mappertest.cc"

int main( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  std::stringstream gridFile;
  gridFile << GRIDDIM<<"dgrid.dgf";

  Dune::Fem::LagrangeBasis_Test(gridFile.str()).run();
  Dune::Fem::LagrangeMapper_Test< Dune::GridSelector::GridType >( gridFile.str() ).run();
  return 0;
}
