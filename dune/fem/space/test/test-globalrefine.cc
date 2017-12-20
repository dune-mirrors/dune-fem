#include <config.h>

#include <iostream>
#include <string>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/exactsolution.hh>

template< class GridPart >
void testGlobalRefine ( const std::string &dgfName )
{
  // create grid
  Dune::GridPtr< typename GridPart::GridType > grid( dgfName );
  grid->loadBalance();

  // create grid part
  GridPart gridPart( *grid );

  // create discrete function space
  typedef Dune::Fem::FunctionSpace< typename GridPart::ctype, double, GridPart::dimensionworld, 1 > FunctionSpaceType;
  Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPart, 2 > space( gridPart );

  // create grid function
  Dune::Fem::ExactSolution< FunctionSpaceType > function;
  Dune::Fem::GridFunctionAdapter< decltype( function ), GridPart > gridFunction( "grid function", function, gridPart, 4 );

  // create exact solution
  Dune::Fem::AdaptiveDiscreteFunction< decltype( space ) > uh( "uh", space );
  interpolate( gridFunction, uh );

  // refine grid globally and reinterpolate
  Dune::Fem::GlobalRefine::apply( *grid, 1 );
  interpolate( gridFunction, uh );
}


int main ( int argc, char **argv )
try
{
  typedef Dune::GridSelector::GridType GridType;

  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameters
  Dune::Fem::Parameter::append( argc, argv );

  const std::string dgfName = Dune::Fem::Parameter::getValue< std::string >( "dgf", std::to_string( GridType::dimension ) + "dgrid.dgf" );
  Dune::Fem::Parameter::appendDGF( dgfName );

  // test AdaptiveLeafGridPart
  std::cout << ">>> Testing AdaptiveLeafGridPart..." << std::endl;
  testGlobalRefine< Dune::Fem::AdaptiveLeafGridPart< GridType > >( dgfName );

  // test LeafGridPart
  std::cout << ">>> Testing LeafGridPart..." << std::endl;
  testGlobalRefine< Dune::Fem::LeafGridPart< GridType > >( dgfName );

  return 0;
}
catch( Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
