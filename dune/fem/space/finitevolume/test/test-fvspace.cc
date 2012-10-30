#include <config.h>

// C++ includes
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/fvspace.hh>

#include "../../../test/exactsolution.hh"
#include "../../../test/testgrid.hh"


// range dimension
static const int dimRange = DIMRANGE;

// grid type
typedef Dune::GridSelector::GridType GridType;



// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  { }

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  { }

  std::string prefix () const
  {
    std::stringstream s;
    s << "projection-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};



// algorithm
// ---------

template< class GridType >
double algorithm ( GridType &grid, const int step )
{
  // create grid part
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  // function space type
  typedef typename GridPartType::ctype DomainFieldType;
  static const int dimDomain = GridPartType::dimensionworld;
  typedef Dune::Fem::FunctionSpace< DomainFieldType, double, dimDomain, dimRange > FunctionSpaceType;

  // create exact solution
  typedef Dune::Fem::ExactSolution< FunctionSpaceType > ExactSolutionType;
  ExactSolutionType exactSolution;
  typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", exactSolution, gridPart, 5 );

  // create discrete function space
  typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  // create discrete function
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  // perform the L2Projection
  Dune::Fem::L2Projection< GridExactSolutionType, DiscreteFunctionType > dgl2;
  dgl2( gridExactSolution, solution );

  // prepare output
  typedef Dune::tuple< const DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  IOTupleType ioTuple( &solution, &gridExactSolution );

  // data output
  typedef Dune::Fem::DataOutput< GridType, IOTupleType > DataOutputType;
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );
  dataOutput.write();

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  return l2norm.distance( gridExactSolution, solution );
}



int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " [initial refinement level] [steps]" << std::endl;
    return 0;
  }

  // initial refinement level 
  int refCount = 1;
  if( argc > 1 )
    refCount = std::max( refCount, std::atoi( argv[ 1 ] ) );

  // number of EOC steps
  int steps = 1;
  if( argc > 2 )
    steps = std::max( steps, std::atoi( argv[ 2 ] ) );

  // create grid
  GridType &grid = Dune::Fem::TestGrid::grid();
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf(); 
  Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

  // compute DG L2-projection
  double oldError = algorithm( grid, 0 );
  std::cout << "L2 error[ 0 ] = " << oldError << std::endl;
  for( int step = 1; step < steps ; ++step )
  {
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    const double newError = algorithm( grid, step );
    const double eoc = log( oldError / newError ) / M_LN2;
    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      std::cout << "L2 error[ " << step << " ] = " << newError << std::endl;
      std::cout << "Eoc[ " << step << " ] = " << eoc << std::endl;
    }
    oldError = newError;
  }
  
  return 0;
}
catch( Dune::Exception e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
