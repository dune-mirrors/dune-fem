#include <config.h>

#include <sstream>
#include <string>
#include <utility>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/function/hierarchical.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/linear/hierarchical.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin/space.hh>
#include <dune/fem/space/discontinuousgalerkin/tuple.hh>

#include <dune/fem/test/exactsolution.hh>
#include <dune/fem/test/massoperator.hh>


// dgfUnitCube
// -----------

inline static std::string dgfUnitCube ( int dimWorld, int cells )
{
  std::string dgf = "DGF\nINTERVAL\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += " 0";
  dgf += "\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += " 1";
  dgf += "\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += (" " + std::to_string( cells ));
  dgf += "\n#\n";
  return dgf;
}



// FunctionSpace
// -------------

template< class GridPart, int dimRange >
using FunctionSpace = Dune::Fem::FunctionSpace< typename GridPart::ctype, double, GridPart::dimensionworld, dimRange >;



// TaylorHoodDGSpace
// -----------------

template< class GridPart >
using VelocityDGSpace = Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpace< GridPart, GridPart::dimensionworld >, GridPart, 2 >;

template< class GridPart >
using PressureDGSpace = Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpace< GridPart, 1 >, GridPart, 1 >;

template< class GridPart >
using TaylorHoodDGSpace = Dune::Fem::TupleDiscontinuousGalerkinSpace< VelocityDGSpace< GridPart >, PressureDGSpace< GridPart > >;



// main
// ----

int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // construct unit cube
  typedef typename Dune::GridSelector::GridType GridType;
  std::istringstream dgf( dgfUnitCube( GridType::dimensionworld, 16 ) );
  Dune::GridPtr< GridType > grid( dgf );

  // create leaf grid part
  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  // construct Taylor-Hood DG space
  TaylorHoodDGSpace< GridPartType > dfSpace( gridPart );

  // construct exact solution
  Dune::Fem::ExactSolution< FunctionSpace< GridPartType, GridPartType::dimensionworld+1 > > uExact;
  const auto uGridExact = gridFunctionAdapter( "exact solution", uExact, gridPart, 3 );

  // construct discrete function for solution
  typedef Dune::Fem::HierarchicalDiscreteFunction< decltype( dfSpace ) > DiscreteFunctionType;
  DiscreteFunctionType u( "solution", dfSpace );
  u.clear();

  // assemble mass operator
  typedef Dune::Fem::HierarchicalLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  MassOperator< DiscreteFunctionType, LinearOperatorType > massOp( dfSpace );

  // assemble right hand side
  DiscreteFunctionType rhs( "right hand side", dfSpace );
  massOp.assembleRHS( uGridExact, rhs );

  // solve the system using dune-istl
  typedef std::decay_t< decltype( u.dofVector().array() ) > X;
  X &x = u.dofVector().array();
  X &b = rhs.dofVector().array();

  typedef std::decay_t< decltype( massOp.exportMatrix() ) > MatrixType;
  Dune::MatrixAdapter< MatrixType, X, X > A( massOp.exportMatrix() );
  Dune::SeqJac< MatrixType, X, X, 2 > P( massOp.exportMatrix(), 1, 1 );

  Dune::CGSolver< X > solver( A, P, 1e-10, 1000, 0 );
  Dune::InverseOperatorResult result;
  solver.apply( x, b, result );

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  std::cout << "L2-Error (ISTL-CG): " << l2norm.distance( uGridExact, u ) << std::endl;

  // solve again using dune-fem CG inverse operator
  u.clear();
  massOp.assembleRHS( uGridExact, rhs );
  Dune::Fem::Solver::CGInverseOperator< DiscreteFunctionType > invOp( massOp, 1e-10, 1e-10, 1000, false );
  invOp( rhs, u );
  std::cout << "L2-Error (FEM-CG):  " << l2norm.distance( uGridExact, u ) << std::endl;

  return 0;
}
