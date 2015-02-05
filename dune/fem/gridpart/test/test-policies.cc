#include <config.h>

#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/gridpart/levelgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>

#include "../../test/testgrid.hh"

// isFake
// ------

template< class GridPart >
std::true_type __isFake ( const Dune::Fem::GridPart2GridViewImpl< GridPart > & );

std::false_type __isFake ( ... );

template< class GridView >
struct isFake
{
  static const bool v = decltype( __isFake( std::declval< typename GridView::GridViewImp >() ) )::value;
};



// type
// ----

template< class GridView >
std::string type ( const GridView &gridPart )
{
  return isFake< GridView >::v ? "fake" : "real";
}



// gridView
// --------

template< class GridPart >
typename GridPart::GridViewType gridView ( const GridPart &gridPart )
{
  return typename GridPart::GridViewType( gridPart );
}



int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  typedef Dune::YaspGrid< 2 > GridType;
  GridType &grid = Dune::Fem::TestGrid::grid();

  // create simple (level) grid part which has underlying grid view
  typedef Dune::Fem::LevelGridPart< GridType > HostGridPartType;
  HostGridPartType hostGridPart( grid, 0 );
  std::cout << "LevelGridPart has " << type( gridView( hostGridPart ) )
            << " grid view (should be \"real\")."  << std::endl;

  // create meta grid part which needs to fake grid view
  typedef Dune::Fem::IdGridPart< HostGridPartType > GridPartType;
  GridPartType gridPart( hostGridPart );
  std::cout << "IdGridPart has " << type( gridView( gridPart ) )
            << " grid view (should be \"fake\")." << std::endl;

  // range-based for loop using dune-grid facilities
  const GridPartType::IndexSetType &indexSet = gridPart.indexSet();
  for( const auto &entity : Dune::elements( gridView( gridPart ), Dune::Partitions::interiorBorder ) )
  {
    if( !indexSet.contains( entity ) )
      DUNE_THROW( Dune::InvalidStateException, "" );
  }

  return 0;
}
catch( Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
}
catch( std::exception &exception )
{
  std::cerr << exception.what() << std::endl;
}
catch( ... )
{
  std::cerr << "Unknown exception thrown" << std::endl;
}
