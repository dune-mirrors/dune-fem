#ifndef  DUNE_FEM_P12DSPACE_TEST_MAIN_CC__
#define  DUNE_FEM_P12DSPACE_TEST_MAIN_CC__

#include  <config.h>

#include  <dune/fem/misc/suite.hh>

#include  "test_l2projection.hh"
#include  "test_poisson.hh"

using namespace Dune;

int main(int argc, char *argv[])
{
  MPIManager :: initialize( argc, argv );

  std::string gridFile( "2dgrid.dgf" );

  Dune::Suite suite("Generic Basefunction tests");
  suite.addTest( new Dune::L2Projection_Test< GridType >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< GridType >(gridFile) );

  suite.run();
  suite.report();
}

#endif  /*DUNE_FEM_P12DSPACE_TEST_MAIN_CC__*/
