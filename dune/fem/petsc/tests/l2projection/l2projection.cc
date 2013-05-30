// include configurations options

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_PETSC

// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>


// Includes from DUNE-FEM
// ----------------------

// include Lagrange discrete function space
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/petscsolver.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/femeoc.hh>


#include "massoperator.hh"
#include "base.hh"



// Global Type Definitions
// -----------------------

#if defined POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif



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
  typedef Dune::GridSelector :: GridType GridType;
  //typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
  typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > SpaceType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, polOrder > DiscreteSpaceType;
#ifdef PETSCLINEAROPERATOR
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
#else 
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
#endif
  typedef Dune::Fem::L2Norm< GridPartType > L2NormType;

  typedef Function< SpaceType > FunctionType;

  explicit Algorithm ( GridType &grid );
  void operator() ( DiscreteFunctionType &u );
  void finalize ( DiscreteFunctionType &u );
  DiscreteSpaceType &space ();

private:
  GridPartType gridPart_;
  DiscreteSpaceType dfSpace_;
  FunctionType function_;
  int eocId_;
};

inline Algorithm::Algorithm ( GridType &grid )
: gridPart_( grid ),
  dfSpace_( gridPart_, Dune::All_All_Interface )
{
  std::vector< std::string > eocHeaders( 2 );
  eocHeaders[ 0 ] = "$L^2$-Error";
  eocHeaders[ 1 ] = "$H^1$-Error";
  eocId_ = Dune::Fem::FemEoc::addEntry( eocHeaders );
}


inline void Algorithm::operator() ( DiscreteFunctionType &solution )
{
  // get operator
  typedef MassOperator< DiscreteFunctionType > MassOperatorType;
  MassOperatorType massOperator( dfSpace_ );

  // assemble RHS
  DiscreteFunctionType rhs( "rhs", dfSpace_ );
  massOperator.assembleRHS( function_, rhs );

  // get solver
  typedef Dune::Fem::Operator< DiscreteFunctionType > InvOperatorType; 
  InvOperatorType* op = 0;

  unsigned long maxIter = dfSpace_.size();
  maxIter = dfSpace_.grid().comm().sum( maxIter );

#ifdef USE_PETSCLINEAROPERATOR 
  if( Dune::Fem::Parameter::getValue<bool>("usepetsc", true ) )
  {
    typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, MassOperatorType::LinearOperatorType > InverseOperator;  
    op = new InverseOperator( massOperator.systemMatrix(), 1e-10, 1e-10, maxIter );
  }
  else 
#endif
  {
    typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > InverseOperator;  
    op = new InverseOperator( massOperator.systemMatrix(), 1e-10, 1e-10, maxIter );
  }

  // clear result
  solution.clear();

  // apply solver
  InvOperatorType& inverseOperator = *op;
  inverseOperator( rhs, solution );

  Dune::Fem::VTKIO< GridPartType > vtkio( dfSpace_.gridPart() );
  vtkio.addVertexData( solution );
  vtkio.write("dump");

  delete op;
}


inline void Algorithm::finalize ( DiscreteFunctionType &solution )
{
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridFunctionType;

  const int order = DiscreteSpaceType::polynomialOrder+1;
  GridFunctionType gridFunction( "exact solution", function_, gridPart_, order );

  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart_ );
  Dune::Fem::H1Norm< GridPartType > h1norm( gridPart_ );

  std::vector< double > errors( 2 );
  errors[ 0 ] = l2norm.distance( gridFunction, solution );
  errors[ 1 ] = h1norm.distance( gridFunction, solution );

  Dune::Fem::FemEoc::setErrors( eocId_, errors );
}


inline Algorithm::DiscreteSpaceType &Algorithm::space ()
{
  return dfSpace_;
}
#endif


// Main Program
// ------------

int main ( int argc, char **argv )
try
{
#if HAVE_PETSC
  typedef Dune::GridSelector::GridType GridType;

  // initialize MPI manager and PETSc
  Dune::Fem::MPIManager::initialize( argc, argv );

  // add command line parameters to global parameter table
  Dune::Fem::Parameter::append( argc, argv );
  // append parameters from the parameter file
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  // create the grid
  Dune::GridPtr< GridType > gridptr = initialize< GridType >( std::string( "L2 projection" ) );

  Algorithm algorithm( *gridptr );
  compute( algorithm );

  Dune::Fem::Parameter::write( "parameter.log" );
#endif

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
