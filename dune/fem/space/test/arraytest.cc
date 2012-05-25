#include <config.h>

#include <iostream>
#include <dune/common/stdstreams.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

using namespace Dune;

#include <dune/fem/space/common/arrays.hh>

template <class Vector> 
void checkVector( Vector& vector , const size_t rows, const size_t cols ) 
{
  std::cout << "Start Vector check: " << std::endl;
  vector.resize( rows );
  std::cout << "vector: capacity " << vector.capacity() << "  size " << vector.size() << std::endl;
  for( size_t i=0; i<rows; ++i ) 
  {
    vector[ i ].resize( cols );
    std::cout << "vector["<<i<<"]: capacity " << vector[ i ].capacity() << "  size " << vector[ i ].size() << std::endl;
    typedef typename Vector :: value_type :: value_type value_type;
    for( size_t j=0; j<cols; ++j ) 
    {
      vector[ i ][ j ] = value_type( j );
    }
  }
  std::cout << "End Vector check!" << std::endl << std::endl;
}

//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main( int argc, char *argv[] )
try {
  typedef FieldVector< double , 1 > ValueType;

  MutableArray< MutableArray< ValueType > > vector;
  vector.setMemoryFactor( 1.1 );

  checkVector( vector, 10, 4 );
  checkVector( vector, 7, 6 );
  checkVector( vector, 20, 9 );
  checkVector( vector, 2, 5 );
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}

