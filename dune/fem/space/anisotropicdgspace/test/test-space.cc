#include <config.h>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/anisotropicdgspace.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include "../../../test/testgrid.hh"


// get parameters
static const int maxOrder = POLORDER;
static const int dimRange = DIMRANGE;

// grid type
typedef Dune::GridSelector::GridType GridType;
// grid part type
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

// domain field type
typedef GridPartType::ctype DomainFieldType;
// range field type
typedef double RangeFieldType;
// world dimension
static const int dimDomain = GridPartType::dimensionworld;
typedef Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

// anisotropic DG space
typedef Dune::Fem::AnisotropicDGSpace< FunctionSpaceType, GridPartType, maxOrder > DiscreteFunctionSpaceType;

// discrete function space
// typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;


int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  GridType &grid = Dune::TestGrid::grid();
  const int step = Dune::TestGrid::refineStepsForHalf();
  grid.globalRefine( 2*step );

  // create grid part
  GridPartType gridPart( grid );

  // create discrete function space
  typedef DiscreteFunctionSpaceType::MultiIndexType MultiIndexType;
  MultiIndexType multiIndex( maxOrder );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart, multiIndex );

  // create discrete function
  // DiscreteFunctionType uh( "uh", discreteFunctionSpace );

  return 0;
}
catch( Dune::Exception e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
