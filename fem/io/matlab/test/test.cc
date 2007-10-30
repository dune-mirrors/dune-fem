#include <config.h>

#include <iostream>
#include <dune/grid/common/gridpart.hh>

#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/io/matlab/matlabstream.hh>

#include "testgrid.hh"

using namespace Dune;

typedef LeafGridPart< GridType > GridPartType; 

int main()
{
  DynamicVector< double > v( 3 );
  for( unsigned int i = 0; i < 3; ++i )
    v[ i ] = i * (i+1);
  std :: cout << v << std :: endl;

  MatlabOutStream s( "test.vector" );
  s << v;

  DenseMatrix<double> m(4,3);
  for ( unsigned int i =0; i< 4;++i)
      {
	for(unsigned int j =0 ; j < 3; ++j){
          m[i][j]=j*i+j;
       }
      }
  m.print(std::cout);
  s << m;

  GridType &grid = TestGrid :: grid();
  GridPartType gridPart( grid );
  MatlabOutStream sg( "test.grid" );
  sg << gridPart;
}
