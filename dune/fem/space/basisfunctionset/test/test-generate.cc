#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/space/basisfunctionset/codegen.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/simple.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>

#include "checkbasisfunctionset.hh"
#include <dune/fem/test/testgrid.hh>

template< class GridPartType, int polorder >
void generateCode ( GridPartType &gridPart )
{
  static const int dimDomain = GridPartType::dimensionworld;

  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, dimDomain, 1 > ScalarFunctionSpaceType;

  std::vector< int > elemQuadOrds = {{ 3, 4, 5 }};
  std::vector< int > faceQuadOrds = {{ 3, 4 }};
  const std::string path = ".";

  {
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< ScalarFunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
    DiscreteFunctionSpaceType space( gridPart );
    Dune::Fem::generateCode( space, elemQuadOrds, faceQuadOrds, path );
  }
#if HAVE_DUNE_LOCALFUNCTIONS
  {
    typedef Dune::Fem::LagrangeSpace< ScalarFunctionSpaceType, GridPartType, Dune::GaussLobattoPointSet > DiscreteFunctionSpaceType;
    DiscreteFunctionSpaceType space( gridPart, polorder );
    Dune::Fem::generateCode( space, elemQuadOrds, faceQuadOrds, path );
  }
#endif
  {
    typedef Dune::Fem::DiscontinuousGalerkinSpace< ScalarFunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
    DiscreteFunctionSpaceType space( gridPart );
    Dune::Fem::generateCode( space, elemQuadOrds, faceQuadOrds, path );
  }
  {
    typedef Dune::Fem::LegendreDiscontinuousGalerkinSpace< ScalarFunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
    DiscreteFunctionSpaceType space( gridPart );
    Dune::Fem::generateCode( space, elemQuadOrds, faceQuadOrds, path );
  }
}

int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( argc >= 2 ? argv[ 1 ] : "parameter" );

  typedef Dune::GridSelector::GridType GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();

  grid.globalRefine( 1 );

  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  if( gridPart.begin< 0 >() == gridPart.end< 0 >() )
    return 1;

  generateCode< GridPartType, 2 >( gridPart );
  generateCode< GridPartType, 3 >( gridPart );
  generateCode< GridPartType, 4 >( gridPart );

  return 0;
}
