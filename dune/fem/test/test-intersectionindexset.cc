#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

template< class GridPart >
void checkIntersectionIndexSet ( const GridPart &gridPart )
{
  const auto &indexSet = gridPart.indexSet();

  std::size_t count = 0;
  for( const auto &element : elements( gridPart, Dune::Partitions::interiorBorder ) )
  {
    for( const auto &intersection : intersections( gridPart, element ) )
    {
      indexSet.index( intersection );
      ++count;
      if( !intersection.neighbor() )
        ++count;
    }
  }
  count /= 2;

  const std::size_t size = indexSet.size( GridPart::dimension+1 );
  if( count != size )
    DUNE_THROW( Dune::Exception, "IntersectionIndexSet reports wrong size: " << size << " (should be " << count << ")" );
}

const char *unitSquare
   = "DGF\n\n"
     "INTERVAL\n"
     "0 0\n"
     "1 1\n"
     "10 10\n"
     "#\n";

// main program
int main(int argc, char ** argv)
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > GridType;
  std::istringstream dgf( unitSquare );
  Dune::GridPtr< GridType > grid( dgf );

  Dune::Fem::IntersectionAdaptiveLeafGridPart< GridType > gridPart( *grid );
  checkIntersectionIndexSet( gridPart );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
