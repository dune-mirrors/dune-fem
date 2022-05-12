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

#if HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscinverseoperators.hh>
#elif HAVE_DUNE_ISTL && defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
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
#include <dune/fem/solver/krylovinverseoperators.hh>
#endif
#endif

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/common/bindguard.hh>

#include <dune/fem/test/testgrid.hh>

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
typedef Dune::Fem::ISTLInverseOperator< DiscreteFunctionType, Dune::Fem::SolverParameter::cg > InverseOperatorType;
#elif HAVE_EIGEN && defined USE_EIGEN
typedef Dune::Fem::EigenVector< double > DofVectorType;
typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteSpaceType, DofVectorType > > DiscreteFunctionType;
typedef Dune::Fem::EigenLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::EigenInverseOperator< DiscreteFunctionType > InverseOperatorType;
#else
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
#if HAVE_SUITESPARSE_UMFPACK && defined USE_UMFPACK
typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType > InverseOperatorType;
#else
typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > InverseOperatorType;
#endif
#endif


// Algorithm
// ---------

struct Algorithm
{
  typedef Dune::FieldVector< double, 1 > ErrorType;

  explicit Algorithm ( GridType &grid );
  ErrorType operator() ( int step );
  ErrorType finalize ( DiscreteFunctionType &u );
  DiscreteSpaceType &space ();
  void nextMesh ();

private:
  GridType &grid_;
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
  solution.clear();
  DiscreteFunctionType rhs( "rhs", space );
  rhs.clear();

  // get operator
  LinearOperatorType identityOperator( "identity", space, space );
  identityOperator.reserve( Dune::Fem::DiagonalStencil<DiscreteSpaceType,DiscreteSpaceType>( space, space ) );
  identityOperator.clear();

  // first set all zeros to tests whether switching between SetContrib and
  // AddContrib works
  {
    Dune::Fem::SetLocalContribution< LinearOperatorType > localMatrix( identityOperator );
    typedef Dune::Fem::CachingQuadrature< typename DiscreteSpaceType::GridPartType, 0 > QuadratureType;
    Dune::DynamicVector< typename SpaceType::RangeType > values;
    for( const auto &entity : space )
    {
      auto guard = bindGuard( localMatrix, entity, entity );
      const auto &basis = localMatrix.domainBasisFunctionSet();
      const unsigned int numBasisFunctions = basis.size();
      values.resize( numBasisFunctions, 0. );
      for( const auto qp : QuadratureType( entity, localMatrix.order() ) )
        localMatrix.axpy( qp, values );
    }
  }
  {
    Dune::Fem::AddLocalContribution< LinearOperatorType > localMatrix( identityOperator );
    typedef Dune::Fem::CachingQuadrature< typename DiscreteSpaceType::GridPartType, 0 > QuadratureType;
    Dune::DynamicVector< typename SpaceType::RangeType > values;
    for( const auto &entity : space )
    {
      auto guard = bindGuard( localMatrix, entity, entity );
      const auto &basis = localMatrix.domainBasisFunctionSet();
      const unsigned int numBasisFunctions = basis.size();
      values.resize( numBasisFunctions, 1. );
      for( const auto qp : QuadratureType( entity, localMatrix.order() ) )
        localMatrix.axpy( qp, values );
    }
  }
  // this should not be necessary since it should happen in the local
  // contribution destructor
  identityOperator.finalize();

  std::vector<unsigned int> rows( {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19} );
  std::vector<unsigned int> auxRows;
  identityOperator.setUnitRows( rows, auxRows );
  auto rIt = rhs.dbegin();
  for ( unsigned int r=0; r<rows.size(); ++r)
  {
    *rIt = r;
    ++rIt;
  }

  // identityOperator.finalize(); // finished setting up matrix
  identityOperator( rhs, solution );

  int r = 0;
  for ( auto&& s : dofs(solution) ) // note that this doesn't work for petsc
  {
    if ( find (rows.begin(), rows.end(), r) == rows.end() )
      s = 0;
    ++r;
  }

  solution -= rhs;

  return finalize( solution );
}


inline Algorithm::ErrorType Algorithm::finalize ( DiscreteFunctionType &solution )
{
  ErrorType error( solution.scalarProductDofs( solution ) );
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

  const int nrSteps = 1;

  std::cout<< "testing with polorder "<< polOrder <<std::endl;
  Algorithm algorithm( grid );
  Algorithm::ErrorType *error;
  error = new Algorithm::ErrorType[ nrSteps ];
  for( int step = 0; step<nrSteps; ++step )
  {
    error[ step ] = algorithm( step );
    algorithm.nextMesh();
    std::cout << "result: " << step << " " << error[step] << std::endl;
  }

  for( int step = 0; step <nrSteps; ++step )
  {
    if( std::abs( error[step] )  > 1e-10 )
    {
      DUNE_THROW(Dune::InvalidStateException,"EOC check of solving laplace matrix system failed");
    }
  }

  delete [] error;
  Dune::Fem::Parameter::write( "parameter.log" );

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
