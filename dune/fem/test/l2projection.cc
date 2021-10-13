// include configurations options

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>

//#define USE_OLD_RANNACHERTUREK_SPACE

// Includes from DUNE-FEM
// ----------------------

// include Lagrange discrete function space
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/function/localfunction/bindable.hh>

#if HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscinverseoperators.hh>
#elif HAVE_DUNE_ISTL && defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#if HAVE_SUITESPARSE_UMFPACK && defined USE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#else
#include <dune/fem/solver/krylovinverseoperators.hh>
#endif
#endif

#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>

#include "massoperator.hh"
#include "testgrid.hh"

// Global Type Definitions
// -----------------------

typedef Dune::GridSelector :: GridType GridType;

typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
typedef Dune::Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype, GridType::dimensionworld, 1 > SpaceType;
typedef Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< SpaceType, GridPartType > DiscreteSpaceType;
//typedef Dune::Fem::LagrangeSpace< SpaceType, GridPartType > DiscreteSpaceType;

#if HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
typedef Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType > InverseOperatorType;
#elif HAVE_DUNE_ISTL && defined USE_BLOCKVECTORFUNCTION
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::ISTLInverseOperator< DiscreteFunctionType, Dune::Fem::SolverParameter::cg > InverseOperatorType;
#else
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
#if HAVE_SUITESPARSE_UMFPACK && defined USE_UMFPACK
typedef Dune::Fem::UMFPACKInverseOperator< DiscreteFunctionType > InverseOperatorType;
#else
typedef Dune::Fem::KrylovInverseOperator< DiscreteFunctionType > InverseOperatorType;
#endif
#endif

typedef MassOperator< DiscreteFunctionType, LinearOperatorType > MassOperatorType;

// Function to Project
// -------------------
template <class GridPart, class RangeType>
struct Function : public Dune::Fem::BindableGridFunction< GridPart, RangeType >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, RangeType > Base;
  using Base::Base;

  template <class Point>
  void evaluate ( const Point &p, RangeType &value ) const
  {
    auto x = Base::global(p);
    value = 1.;
    for( int k = 0; k < x.dimension; ++k )
      value *= sin( M_PI * x[k] );
  }
  template <class Point>
  void jacobian ( const Point &p, typename Base::JacobianRangeType &jacobian ) const
  {
    auto x = Base::global(p);
    for( int j = 0; j < x.dimension; ++j )
    {
      // jacobian has only one row, calc j-th column
      jacobian[0][j] = M_PI;
      for( int k = 0; k < x.dimension; ++k )
        jacobian[0][j] *= (j == k ? cos( M_PI*x[k] ) : sin( M_PI*x[k] ));
    }
  }
  template <class Point>
  void hessian ( const Point &p, typename Base::HessianRangeType &h ) const
  {
  }

  unsigned int order() const { return 4; }
  std::string name() const { return "MyFunction"; } // needed for output
};

// Algorithm
// ---------

struct Algorithm
{
  typedef Dune::FieldVector< double, 2 > ErrorType;
  typedef Function< GridPartType, SpaceType::RangeType > FunctionType;
  typedef Dune::FemPy::VirtualizedGridFunction< GridPartType, typename SpaceType::RangeType > VGridFunctionType;

  explicit Algorithm ( GridType &grid );
  ErrorType operator() ( int step, int polOrder );
  ErrorType finalize ( DiscreteFunctionType &u );
  DiscreteSpaceType &space ();
  void nextMesh ();
  void resetMesh (const int steps);

private:
  GridPartType gridPart_;
  FunctionType function_;
};

inline Algorithm::Algorithm ( GridType &grid )
: gridPart_( grid ), function_(gridPart_)
{
}

inline void Algorithm :: nextMesh ()
{
  Dune::Fem::GlobalRefine::apply(gridPart_.grid(), Dune::DGFGridInfo< GridType >::refineStepsForHalf() );
}

inline void Algorithm :: resetMesh (const int steps)
{
  gridPart_.grid().globalRefine( -steps * Dune::DGFGridInfo< GridType >::refineStepsForHalf() );
}

inline Algorithm::ErrorType Algorithm::operator() ( int step, int polOrder )
{
  DiscreteSpaceType space( gridPart_, polOrder );
  DiscreteFunctionType solution( "solution", space );

  // get operator
  MassOperatorType massOperator( space );

  // assemble RHS
  DiscreteFunctionType rhs( "rhs", space );

  // create virtualized grid function to test this class
  VGridFunctionType vF( function_ );

  massOperator.assembleRHS( vF, rhs );

  unsigned long maxIter = space.size();
  maxIter = space.gridPart().comm().sum( maxIter );

  // clear result
  solution.clear();

  // apply solver
  InverseOperatorType inverseOperator;
  inverseOperator.bind( massOperator );

  inverseOperator.parameter().setTolerance( 1e-10 );
  inverseOperator.parameter().setMaxIterations( maxIter );

  inverseOperator( rhs, solution );
  return finalize( solution );
}


inline Algorithm::ErrorType Algorithm::finalize ( DiscreteFunctionType &solution )
{
  const GridPartType &gridPart = solution.space().gridPart();
  ErrorType error;

  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );

  // create virtualized grid function to test this class
  VGridFunctionType vF( function_ );

  DiscreteFunctionType rhs( solution );

  Dune::Fem::interpolate( vF, rhs );

  {
    // create virtualized grid function to test this class
    // this is only for testing, error needs to be recomputed
    VGridFunctionType vF2( rhs );
    error[ 0 ] = l2norm.distance( vF2, solution );
    error[ 1 ] = h1norm.distance( vF2, solution );
  }

  error[ 0 ] = l2norm.distance( vF, solution );
  error[ 1 ] = h1norm.distance( vF, solution );
  std::cout << error <<std::endl;
  return error;
}

void run( GridType& grid, const int polOrder )
{
  const int nrSteps = 4;

  std::cout<< "testing with polorder "<< polOrder <<std::endl;
  Algorithm algorithm( grid );
  std::vector< typename Algorithm::ErrorType > error( nrSteps );
  for( int step = 0; step<nrSteps; ++step )
  {
    error[ step ] = algorithm( step, polOrder );
    algorithm.nextMesh();
  }

  for( int step = 1; step <nrSteps; ++step )
  {
    double l2eoc = log( error[ step ][ 0 ] / error[ step -1 ][ 0 ] ) / log( 0.5 );
    double h1eoc = log( error[ step ][ 1 ] / error[ step -1 ][ 1 ] ) / log( 0.5 );

    if( std::abs( l2eoc -1 - polOrder )  > 0.2 )
    {
      DUNE_THROW(Dune::InvalidStateException,"EOC check of solving mass matrix system failed L2eoc " << l2eoc << " " << polOrder);
    }

    if( std::abs( h1eoc - polOrder )  > 0.2 )
    {
      //note: This will fail with Yaspgrid, bug in Geometry JacobianInverse
      DUNE_THROW(Dune::InvalidStateException,"EOC check of solving mass matrix system failed H1eoc " << h1eoc << " " << polOrder);
    }
  }

  algorithm.resetMesh(nrSteps);
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
  for( int p=1; p<2; ++p )
    run( grid, p );

  Dune::Fem::Parameter::write( "parameter.log" );

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
