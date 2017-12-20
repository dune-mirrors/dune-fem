#include <config.h>

#include "checkgridpart.hh"

int main ( int argc, char ** argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // create grid
  auto& grid = Dune::Fem::TestGrid::grid();

  // refine grid
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();
  grid.globalRefine( 2*step );
  grid.loadBalance();

  // create grid part
  typedef Dune::GridSelector::GridType GridType;
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  typedef typename GridPartType :: GridViewType GridViewType;
  //GridPartType gridPart( grid );
  std::unique_ptr< GridPartType > gpPtr;
  std::unique_ptr< GridViewType > gvPtr;

  gpPtr.reset( new GridPartType( grid ) );
  GridPartType gp2( grid );

  gvPtr.reset( new GridViewType( static_cast<GridViewType> ( *gpPtr) ) );

  gpPtr.reset( new GridPartType( grid ) );
  gvPtr.reset( new GridViewType( static_cast<GridViewType> ( *gpPtr) ) );
  gpPtr.reset();
  gvPtr.reset();

  GridPartType& gridPart = gp2;

  // run tests
  std::cout << "Testing entities" << std::endl;
  testGridPart( gridPart );
  std::cout << std::endl;

  std::cout << "Testing subentities" << std::endl;
  testSubEntities< GridType::dimension >( gridPart );
  std::cout << std::endl;

  std::cout << "GridWidth: " << Dune::Fem::GridWidth::calcGridWidth( gridPart ) << std::endl;



  typedef Dune::DefaultFailureHandler FailureHandlerType;
  FailureHandlerType failureHandler;
  std::cout << "Testing entity seeds" << std::endl;
  Dune::Fem::CheckEntitySeed< GridPartType >::check( gridPart );
  std::cout << "Testing geometries" << std::endl;
  Dune::Fem::CheckGeometry< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing index set" << std::endl;
  Dune::Fem::CheckIndexSet< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );
  std::cout << "Testing intersections" << std::endl;
  Dune::Fem::CheckIntersections< GridPartType, FailureHandlerType >::check( gridPart, failureHandler );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
