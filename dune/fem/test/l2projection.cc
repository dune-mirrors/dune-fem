// include configurations options

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>


#if defined HAVE_PETSC

// Includes from DUNE-FEM
// ----------------------

// include Lagrange discrete function space
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/solver/petscsolver.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include "massoperator.hh"
#include "testgrid.hh"


// Global Type Definitions
// -----------------------

#if defined POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif


typedef Dune::GridSelector :: GridType GridType;


// Function to Project
// -------------------

template< class FunctionSpace >
struct Function
{
  typedef FunctionSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  void evaluate ( const DomainType &x, RangeType &value ) const
  {
    value = 1.;
    for( int k = 0; k < FunctionSpace::dimDomain; ++k )
      value *= sin( M_PI * x[k] );
  }

  void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
  {
    for( int j = 0; j < FunctionSpace::dimDomain; ++j )
    {
      // jacobian has only one row, calc j-th column      
      jacobian[0][j] = M_PI;
      for( int k = 0; k < FunctionSpace::dimDomain; ++k )
        jacobian[0][j] *= (j == k ? cos( M_PI*x[k] ) : sin( M_PI*x[k] ));
    }
  }
};


// Algorithm
// ---------

struct Algorithm
{
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > SpaceType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, polOrder > DiscreteSpaceType;
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;

  typedef Dune::FieldVector< double, 2 > ErrorType;

  typedef Function< SpaceType > FunctionType;

  explicit Algorithm ( GridType &grid );
  ErrorType operator() ( int step );
  ErrorType finalize ( DiscreteFunctionType &u );
  DiscreteSpaceType &space ();
  void nextMesh ();

private:
  GridType &grid_;
  FunctionType function_;
};

inline Algorithm::Algorithm ( GridType &grid )
: grid_( grid )
{
}

inline void Algorithm :: nextMesh () 
{
  grid_.globalRefine( Dune::DGFGridInfo< GridType >::refineStepsForHalf() );
}

inline typename Algorithm::ErrorType Algorithm::operator() ( int step )
{
  GridPartType gridPart( grid_ );
  DiscreteSpaceType space( gridPart );
  DiscreteFunctionType solution( "solution", space );
     
  // get operator
  typedef MassOperator< DiscreteFunctionType > MassOperatorType;
  MassOperatorType massOperator( space );

  // assemble RHS
  DiscreteFunctionType rhs( "rhs", space );
  massOperator.assembleRHS( function_, rhs );

  unsigned long maxIter = space.size();
  maxIter = space.gridPart().comm().sum( maxIter );

  typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, MassOperatorType > InverseOperatorType;  

  // clear result
  solution.clear();

  // apply solver
  InverseOperatorType inverseOperator ( massOperator, 1e-10, 1e-10, maxIter );
  inverseOperator( rhs, solution );

  return finalize( solution );
}


inline typename Algorithm::ErrorType Algorithm::finalize ( DiscreteFunctionType &solution )
{
  const GridPartType &gridPart = solution.space().gridPart();
  ErrorType error;
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridFunctionType;

  const int order = DiscreteSpaceType::polynomialOrder+1;
  GridFunctionType gridFunction( "exact solution", function_, gridPart, order );

  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );

  error[ 0 ] = l2norm.distance( gridFunction, solution );
  error[ 1 ] = h1norm.distance( gridFunction, solution );
  return error;
}

#endif // #if defined HAVE_PETSC

// Main Program
// ------------

int main ( int argc, char **argv )
try
{
#if defined HAVE_PETSC
  // initialize MPI manager and PETSc
  Dune::Fem::MPIManager::initialize( argc, argv );

  // add command line parameters to global parameter table
  Dune::Fem::Parameter::append( argc, argv );
  // append parameters from the parameter file
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  GridType &grid = Dune::Fem::TestGrid :: grid();
  const int step = Dune::Fem::TestGrid :: refineStepsForHalf();
  grid.globalRefine( 2*step );

  const int nrSteps = 4;

  std::cout<< "testing with polorder "<< polOrder <<std::endl;
  Algorithm algorithm( grid );
  typename Algorithm::ErrorType *error; 
  error = new  typename Algorithm::ErrorType[ nrSteps ];
  for( int step = 0; step<nrSteps; ++step )
  {
    error[ step ] = algorithm( step );
    algorithm.nextMesh();
  }

  for( int step = 1; step <nrSteps; ++step )
  {
    // !!! note there is a bug within the H1 Norm, the test will fail !!!
    double l2eoc = log( error[ step ][ 0 ] / error[ step -1 ][ 0 ] ) / log( 0.5 );
    double h1eoc = log( error[ step ][ 1 ] / error[ step -1 ][ 1 ] ) / log( 0.5 );

//    std::cout<< "L2 Eoc: " << l2eoc << std::endl;
//    std::cout<< "H1 Eoc: " << h1eoc << std::endl;
    if( std::abs( l2eoc -1 - polOrder )  > 0.2 )
    {
      DUNE_THROW(Dune::InvalidStateException,"EOC check of refinement failed");
    }
    if( std::abs( h1eoc - polOrder )  > 0.2 )
    { 
      //note: This will fail with Yaspgrid, bug in Geometry JacobianInverse
      DUNE_THROW(Dune::InvalidStateException,"EOC check of refinement failed");
    }
  }

  delete error;
  Dune::Fem::Parameter::write( "parameter.log" );

  return 0;

#else

  std::cout<<"Running an empty Petsc test, so its successful ;)"<<std::endl;
  return 0;

#endif // #if defined HAVE_PETSC
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
