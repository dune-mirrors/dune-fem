#include <config.h>

#include <iostream>

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/io/matlab/matlabstream.hh>

#include "testgrid.hh"

using namespace Dune;


int main()
{
  typedef Dune::GridSelector::GridType GridType;
  typedef LeafGridPart< GridType > GridPartType; 
  // test writing a vector and a matrix to 1 file
  Fem::MatlabOutStream s( "test.data" );
  
  Fem :: DynamicVector< double > v( 3 );
  for( unsigned int i = 0; i < 3; ++i )
    v[ i ] = i * (i+1);
  std :: cout << v << std :: endl;
  s << v;

  Fem::DenseMatrix<double> m(4,3);
  for( unsigned int i =0; i< 4; ++i )
  {
    for( unsigned int j = 0; j < 3; ++j )
      m[ i ][ j ]= j * (i+1);
  }
  m.print( std::cout );
  s << m;

  // test writing out a grid
  GridType &grid = TestGrid :: grid();
  GridPartType gridPart( grid );
  Fem::MatlabOutStream sg( "test.grid" );
  sg << gridPart;


  Fem::MatlabInStream in( "test.vector" );
  in >> v;

}
