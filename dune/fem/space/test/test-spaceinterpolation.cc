#include <config.h>

#include <sstream>
#include <string>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/test/exactsolution.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fem/space/brezzidouglasmarini.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/rannacherturek.hh>

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


static const int polOrder = 2;

typedef Dune::GridSelector::GridType GridType;
typedef Dune::Fem::LeafGridPart< GridType > GridPartType;

static const int dimRange = GridPartType::dimensionworld;

typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype,
        typename GridPartType::ctype, GridPartType::dimensionworld, dimRange > FunctionSpaceType;

typedef Dune::Fem::BrezziDouglasMariniSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType  > DiscreteFunctionSpaceType;
//typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Dune::Fem::LegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Dune::Fem::RannacherTurekSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
//typedef Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;

typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;


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

  for( std::size_t s = 0; s <4; ++s )
  {
    DiscreteFunctionSpaceType space( gridPart );
    DiscreteFunctionType u( "solution", space );

    // interpolate a function
    Dune::Fem::ExactSolution< FunctionSpaceType > uExact;
    const auto uGridExact = gridFunctionAdapter( "exact solution", uExact, gridPart, 3 );
    interpolate( uGridExact, u );

    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
    Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );

    std::cout << l2norm.distance( uGridExact, u ) << "\t" << h1norm.distance( uGridExact, u ) << std::endl;
    grid->globalRefine(1);
  }

  return 0;
}
