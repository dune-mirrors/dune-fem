#include <config.h>

// C++ includes
#include <iostream>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/fvspace.hh>

#include "../../../test/exactsolution.hh"
#include "../../../test/testgrid.hh"


static const int dimRange = DIMRANGE;

// grid type
typedef Dune::GridSelector::GridType GridType;

// grid part type
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

// function space type
typedef GridPartType::ctype DomainFieldType;
typedef double RangeFieldType;
static const int dimDomain = GridPartType::dimensionworld;
typedef Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

// finite volume space type
typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;

// discrete function type
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;


int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  GridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();
  grid.globalRefine( 2*step );

  // create grid part
  GridPartType gridPart( grid );

  // create discrete function space
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  // create exact solution
  typedef Dune::Fem::ExactSolution< FunctionSpaceType > ExactSolutionType;
  ExactSolutionType exactSolution;
  // create discrete function
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  // perform the L2Projection
  Dune::Fem::L2Projection< ExactSolutionType, DiscreteFunctionType > dgl2;
  dgl2( exactSolution, solution );

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  RangeFieldType error = l2norm.distance( exactSolution, solution );
  
  return 0;
}
catch( Dune::Exception e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
