// include configurations options

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>


// Includes from DUNE-FEM
// ----------------------

// include Lagrange discrete function space
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>

#if HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#elif HAVE_DUNE_ISTL && defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#elif HAVE_EIGEN && defined USE_EIGEN
#include <dune/fem/storage/eigenvector.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/solver/eigen.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#if HAVE_SUITESPARSE_UMFPACK && defined USE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#else
#include <dune/fem/solver/cginverseoperator.hh>
#endif
#endif

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

typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
typedef Dune::Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype, GridType::dimensionworld, 1 > SpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, polOrder > DiscreteSpaceType;
#if HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
typedef Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType > InverseOperatorType;
#elif HAVE_DUNE_ISTL && defined USE_BLOCKVECTORFUNCTION
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > InverseOperatorType;
#elif HAVE_EIGEN && defined USE_EIGEN
typedef Dune::Fem::EigenVector< double > DofVectorType;
typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteSpaceType, DofVectorType > > DiscreteFunctionType;
typedef Dune::Fem::EigenLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::EigenCGInverseOperator< DiscreteFunctionType > InverseOperatorType;
#else
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
#if HAVE_SUITESPARSE_UMFPACK && defined USE_UMFPACK
typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType > InverseOperatorType;
#else
typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > InverseOperatorType;
#endif
#endif

typedef MassOperator< DiscreteFunctionType, LinearOperatorType > MassOperatorType;

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

inline Algorithm::ErrorType Algorithm::operator() ( int step )
{
  GridPartType gridPart( grid_ );
  DiscreteSpaceType space( gridPart );
  DiscreteFunctionType solution( "solution", space );

  // get operator
  MassOperatorType massOperator( space );

  // assemble RHS
  DiscreteFunctionType rhs( "rhs", space );
  massOperator.assembleRHS( function_, rhs );

  unsigned long maxIter = space.size();
  maxIter = space.gridPart().comm().sum( maxIter );

  // clear result
  solution.clear();

  // apply solver
  InverseOperatorType inverseOperator ( massOperator, 1e-10, 1e-10, maxIter );
  inverseOperator( rhs, solution );

  return finalize( solution );
}


inline Algorithm::ErrorType Algorithm::finalize ( DiscreteFunctionType &solution )
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


// Main Program
// ------------

int main ( int argc, char **argv )
try
{
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
  Algorithm::ErrorType *error;
  error = new  Algorithm::ErrorType[ nrSteps ];
  for( int step = 0; step<nrSteps; ++step )
  {
    error[ step ] = algorithm( step );
    algorithm.nextMesh();
  }

  for( int step = 1; step <nrSteps; ++step )
  {
    double l2eoc = log( error[ step ][ 0 ] / error[ step -1 ][ 0 ] ) / log( 0.5 );
    double h1eoc = log( error[ step ][ 1 ] / error[ step -1 ][ 1 ] ) / log( 0.5 );

//    std::cout<< "L2 Eoc: " << l2eoc << std::endl;
//    std::cout<< "H1 Eoc: " << h1eoc << std::endl;
    if( std::abs( l2eoc -1 - polOrder )  > 0.2 )
    {
      DUNE_THROW(Dune::InvalidStateException,"EOC check of solving mass matrix system failed");
    }
    if( std::abs( h1eoc - polOrder )  > 0.2 )
    {
      //note: This will fail with Yaspgrid, bug in Geometry JacobianInverse
      DUNE_THROW(Dune::InvalidStateException,"EOC check of solving mass matrix system failed");
    }
  }

  delete error;
  Dune::Fem::Parameter::write( "parameter.log" );

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
