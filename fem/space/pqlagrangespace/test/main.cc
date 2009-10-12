#ifndef  DUNE_FEM_P12DSPACE_TEST_MAIN_CC__
#define  DUNE_FEM_P12DSPACE_TEST_MAIN_CC__

#include  <config.h>

#include  <dune/fem/misc/suite.hh>

#if HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include  "test_l2projection.hh"
#include  "test_poisson.hh"

using namespace Dune;

int main(int argc, char *argv[])
{
  MPIManager :: initialize( argc, argv );

  const std::string gridFile  ("2dgrid.dgf");
  const std::string grid3DFile("3dgrid.dgf");

  // YaspGrid quadrilateral grid
  typedef Dune :: YaspGrid< 2 >                                      QuadGridUnitSquare;
  Dune::Suite suite("Generic Basefunction tests");
  suite.addTest( new Dune::L2Projection_Test< QuadGridUnitSquare, 1, true >(gridFile) );
  suite.addTest( new Dune::L2Projection_Test< QuadGridUnitSquare, 2, true >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< QuadGridUnitSquare, 1, true >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< QuadGridUnitSquare, 2, true >(gridFile) );

#ifdef HAVE_ALUGRID
  // Alugrid simplical grid
  typedef Dune :: ALUSimplexGrid< 2, 2 >                             SimplexGridUnitSquare;
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitSquare, 1 >(gridFile) );
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitSquare, 2 >(gridFile) );
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitSquare, 3 >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitSquare, 1 >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitSquare, 2 >(gridFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitSquare, 3 >(gridFile) );

  typedef Dune :: ALUSimplexGrid< 3, 3 >                             SimplexGridUnitCube;
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitCube, 1 >(grid3DFile) );
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitCube, 2 >(grid3DFile) );
  suite.addTest( new Dune::L2Projection_Test< SimplexGridUnitCube, 3 >(grid3DFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitCube, 1 >(grid3DFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitCube, 2 >(grid3DFile) );
  suite.addTest( new Dune::Poisson_Test< SimplexGridUnitCube, 3 >(grid3DFile) );
#endif

  suite.run();
  suite.report();
}

#endif  /*DUNE_FEM_P12DSPACE_TEST_MAIN_CC__*/
