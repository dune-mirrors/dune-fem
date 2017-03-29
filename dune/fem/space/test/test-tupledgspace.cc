#include <config.h>

#include <sstream>
#include <string>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin/space.hh>
#include <dune/fem/space/discontinuousgalerkin/tuple.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/test/exactsolution.hh>


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



// main
// ----

int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // construct unit cube
  typedef typename Dune::GridSelector::GridType GridType;
  std::istringstream dgf( dgfUnitCube( GridType::dimensionworld, 4 ) );
  Dune::GridPtr< GridType > grid( dgf );

  // create leaf grid part
  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  // construct Taylor-Hood DG space
  typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpace< GridPartType, GridPartType::dimensionworld >, GridPartType, 2 > VelocityDiscreteFunctionSpaceType;
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpace< GridPartType, 1 >, GridPartType, 1 > PressureDiscreteFunctionSpaceType;
  typedef Dune::Fem::TupleDiscontinuousGalerkinSpace< VelocityDiscreteFunctionSpaceType, PressureDiscreteFunctionSpaceType > DiscreteFunctionSpaceType;
  DiscreteFunctionSpaceType dfSpace( gridPart );

  // interpolate a function
  Dune::Fem::ExactSolution< FunctionSpace< GridPartType, GridPartType::dimensionworld+1 > > uExact;
  const auto uGridExact = gridFunctionAdapter( "exact solution", uExact, gridPart, 3 );
  Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > u( "solution", dfSpace );
  interpolate( uGridExact, u );

  // output analytical function and interolation to vtk file
  Dune::Fem::VTKIO< GridPartType > vtkIO( gridPart, Dune::VTK::nonconforming );
  vtkIO.addVertexData( uGridExact );
  vtkIO.addVertexData( u );
  vtkIO.write( "test-tupledgspace" );

  return 0;
}
