// include configurations options
#include <config.h>


// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>


// Includes from Core Modules
// --------------------------

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>


// Includes from DUNE-FEM
// ----------------------

// include support for versioning
#include <dune/fem/version.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// include basic grid parts
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include Lagrange discrete function space
#include <dune/fem/space/pqlagrangespace/pk2dlagrangespace.hh>
#include <dune/fem/space/pqlagrangespace/pq22dlagrangespace.hh>

#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/solver/inverseoperators.hh>

#include <dune/fem/misc/l2error.hh>

#include <dune/fem/io/file/vtkio.hh>

#include "massoperator.hh"



// Global Type Definitions
// -----------------------

#if defined POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

// select grid part to use
typedef Dune::AdaptiveLeafGridPart< Dune::GridSelector::GridType > GridPartType;

typedef Dune::FunctionSpace< double, double, Dune::GridSelector::GridType::dimension, 1 > FunctionSpaceType;

// typedef Dune::Pk2DLagrangeSpace< FunctionSpaceType, GridPartType, polOrder >
//   DiscreteFunctionSpaceType;
typedef Dune::PQ22DLagrangeSpace< FunctionSpaceType, GridPartType >
  DiscreteFunctionSpaceType;

typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;



template< class FunctionSpace >
struct Function
{
  typedef typename FunctionSpace::DomainType DomainType;
  typedef typename FunctionSpace::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &value ) const
  {
    value = 1.;
    for( int i = 0; i < FunctionSpace::dimDomain; ++i )
      value *= sin( M_PI * x[ i ] );
  }

  void evaluate ( const DomainType &x, const double &t, RangeType &value ) const
  {
    evaluate( x, value );
  }
};



// Algorithm
// ---------

template< class Function >
double algorithm ( const Function &function, Dune::GridSelector::GridType &grid, int repeat )
{
  typedef Dune::MassOperator< DiscreteFunctionType > MassOperatorType;
  typedef Dune::CGInverseOp< DiscreteFunctionType, MassOperatorType > InverseOperator;

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType dfSpace( gridPart );

  MassOperatorType massOperator( dfSpace );
  DiscreteFunctionType functional( "functional", dfSpace );
  massOperator( function, functional );

  DiscreteFunctionType solution( "solution", dfSpace );
  solution.clear();
  InverseOperator inverseOperator( massOperator, 1e-10, 1e-10, 50000, 0 );
  inverseOperator( functional, solution );

  if( Dune::Parameter::getValue< bool >( "l2projection.vtkout", false ) )
  {
    Dune::SubsamplingVTKIO< GridPartType > vtkWriter( gridPart, 4 );
    vtkWriter.addVertexData( solution );
    std::ostringstream name;
    name << "l2projection-" << repeat;
    std::cerr << "Writing solution to " << name.str() << "..." << std::endl;
    vtkWriter.write( name.str() );
  }
  
  Dune::L2Error< DiscreteFunctionType > l2error;
  return l2error.norm( function, solution );
}



// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  // initialize MPI manager
  Dune::MPIManager::initialize( argc, argv );
  // add command line parameters to global parameter table
  Dune::Parameter::append( argc, argv );
  Dune::Parameter::append( "parameter" );

  // get name of grid file
  std::string gridfile;
  Dune::Parameter::get( "l2projection.gridfile", gridfile );

  const int level = Dune::Parameter::getValue< int >( "l2projection.level" );
  const int repeats = Dune::Parameter::getValue< int >( "l2projection.repeats", 0 );

  // use DGF parser to create grid
  typedef Dune::GridSelector::GridType GridType;
  Dune::GridPtr< GridType > gridptr( gridfile );
  // distribute the grid to all processes
  gridptr->loadBalance();
  gridptr->globalRefine( level );

  Function< FunctionSpaceType > function;

  double error = algorithm( function, *gridptr, 0 );
  std::cout << error << std::endl;
  for( int i = 0; i < repeats; ++i )
  {
    const double prevError = error;
    gridptr->globalRefine( Dune::DGFGridInfo< GridType >::refineStepsForHalf() );
    error = algorithm( function, *gridptr, i+1 );
    const double eoc = log( prevError / error ) / M_LN2;
    std::cout << error << "  " << eoc << std::endl;
  }

  return 0;
}
catch( Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception.what() << std::endl;
  return 1;
}
