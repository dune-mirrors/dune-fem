#include <config.h>

#define WANT_CACHED_COMM_MANAGER 1

#include <cstdlib>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

#include <dune/common/exceptions.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/test/exactsolution.hh>
#include <dune/fem/test/testgrid.hh>

int main ( int argc, char **argv )
try
{
  using Dune::Fem::Parameter;

  typedef Dune::GridSelector::GridType GridType;

  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameters
  Dune::Fem::Parameter::append( argc, argv );

  std::ostringstream s;
  s << GridType::dimension << "dgrid_8.dgf";

  // create grid
  GridType &grid = Dune::Fem::TestGrid::grid( s.str() );

  // initial refinement
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf();
  const int count = Parameter::getValidValue< int >( "startlevel", 1, [] ( int level ) { return (level >= 0); } );
  Dune::Fem::GlobalRefine::apply( grid, count*refineStepsForHalf );

  // create grid part
  //typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  // create grid function
  typedef Dune::Fem::FunctionSpace< GridPartType::ctype, double, GridPartType::dimensionworld, 3 > FunctionSpaceType;
  Dune::Fem::ExactSolution< FunctionSpaceType > function;
  Dune::Fem::GridFunctionAdapter< decltype( function ), GridPartType > gridFunction( "grid function", function, gridPart, 1 );

  // create discrete function space
  typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType space( gridPart );

  // create discrete function
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  DiscreteFunctionType uh( "uh", space );
  Dune::Fem::interpolate( gridFunction, uh );

  // file I/O
  auto tuple = std::make_tuple( &uh );
  Dune::Fem::DataOutput< GridType, decltype( tuple ) > output( grid, tuple );
  output.write();

  // create adaptation manager for grid refinement
  Dune::Fem::RestrictProlongDefaultTuple< DiscreteFunctionType > rp( uh );
  Dune::Fem::AdaptationManager< GridType, decltype( rp ) > adaptManager( grid, rp );

  // generate random seed
  unsigned int seed = 0;
  if( !Parameter::exists( "seed" ) )
  {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    grid.comm().broadcast( &seed, 1, 0 );
  }
  else
    seed = Parameter::getValue< unsigned int >( "seed" );
  if( grid.comm().rank() == 0 )
    std::cerr << "Using seed " << seed << std::endl;

  // create random distribution
  std::default_random_engine engine( seed );
  std::uniform_int_distribution< int > distribution( -(1 << GridType::dimension), 1 );

  const int steps = Parameter::getValidValue< int >( "steps", 3, [] ( int steps ) { return (steps > 0); } );
  for( int step = 0; step < steps; ++step )
  {
    for( const auto &entity : elements( grid.leafGridView() ) )
    {
      auto mark = distribution( engine ) < 0 ? -1 : 1;
      if( entity.level() == 0 )
        mark = std::max( mark, 0 );
      grid.mark( mark, entity );
    }
    adaptManager.adapt();
    output.write();
  }

  return 0;
}
catch( Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
