#include <config.h>

#include <chrono>
#include <iostream>
#include <random>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/solver/timeprovider.hh>

int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // start and end time for time loop and maximum time step
  const double start = 0., end = 1., maxStep = 0.2;

  // test fixed step time provider
  Dune::Fem::FixedStepTimeProvider<> fixedStepTimeProvider(start,maxStep);
  for( ; fixedStepTimeProvider.time() < end; fixedStepTimeProvider.next())
  {
    // print some info
    std::cerr << "FixedStepTimeProvider<>: time step = " << fixedStepTimeProvider.timeStep() << ", "
              << "time = " << fixedStepTimeProvider.time() << ", "
              << "dt = " << fixedStepTimeProvider.deltaT() << std::endl;

  }
  // create random number generator
  const unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine engine( seed );
  std::uniform_real_distribution< double > distribution( start, end );

  // function object returning a random time step estimate
  auto timeStepEstimate = [&distribution, &engine]()
  {
    return distribution( engine );
  };

  // create time provider and provide initial time step estimate
  Dune::Fem::TimeProvider<> timeProvider( start );
  timeProvider.provideTimeStepEstimate( timeStepEstimate() );

  for( timeProvider.init(); timeProvider.time() < end; timeProvider.next() )
  {
    // provide next estimate
    timeProvider.provideTimeStepEstimate( timeStepEstimate() );

    // invalidate time step if estimate was chosen too large and restart
    if( timeProvider.deltaT() > maxStep )
    {
      timeProvider.invalidateTimeStep();
      continue;
    }

    // print some info
    std::cerr << "TimeProvider<>: time step = " << timeProvider.timeStep() << ", "
              << "time = " << timeProvider.time() << ", "
              << "dt = " << timeProvider.deltaT() << std::endl;
  }

  return 0;
}
catch( Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
catch( std::exception &exception )
{
  std::cerr << exception.what() << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Unknown exception thrown" << std::endl;
  return 1;
}
