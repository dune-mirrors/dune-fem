#include <config.h>

#include <iostream>
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/io/matlab/matlabstream.hh>

#include "testgrid.hh"

using namespace Dune;

typedef LeafGridPart< GridType > GridPartType; 

int main()
{
  // test writing a vector and a matrix to 1 file
  MatlabOutStream s( "test.data" );
  
  DynamicVector< double > v( 3 );
  for( unsigned int i = 0; i < 3; ++i )
    v[ i ] = i * (i+1);
  std :: cout << v << std :: endl;
  s << v;

  DenseMatrix<double> m(4,3);
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
  MatlabOutStream sg( "test.grid" );
  sg << gridPart;


  MatlabInStream in( "test.vector" );
  in >> v;

}
