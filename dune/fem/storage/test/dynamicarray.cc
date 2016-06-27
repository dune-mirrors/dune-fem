#include <config.h>

#include <algorithm>
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/storage/dynamicarray.hh>

template <class Vector>
void checkVector( Vector& vector , std::size_t rows, std::size_t cols )
{
  std::cout << "Start Vector check: " << std::endl;
  vector.resize( rows );
  std::cout << "vector: capacity " << vector.capacity() << "  size " << vector.size() << std::endl;
  for( std::size_t i=0; i<rows; ++i )
  {
    vector[ i ].resize( cols );
    std::cout << "vector["<<i<<"]: capacity " << vector[ i ].capacity() << "  size " << vector[ i ].size() << std::endl;
    typedef typename Vector :: value_type :: value_type value_type;
    for( std::size_t j=0; j<cols; ++j )
    {
      vector[ i ][ j ] = value_type( j );
    }
  }
  std::cout << "End Vector check!" << std::endl << std::endl;
}

template <class Vector>
void checkPODVector( Vector& vector , std::size_t rows)
{
  typedef typename Vector :: value_type  value_type;
  if( std::is_same< Dune::Fem::PODArrayAllocator< value_type >, typename Vector::AllocatorType > :: value )
  {
    std::cout << "Using POD allocator for array memory management" << std::endl;
  }

  std::cout << "Start Vector check: " << std::endl;
  vector.resize( rows );
  std::cout << "vector: capacity " << vector.capacity() << "  size " << vector.size() << std::endl;
  std::fill( vector.begin(), vector.end(), value_type( 1 ) );
  std::cout << "End Vector check!" << std::endl << std::endl;
}

int main( int argc, char *argv[] )
try {
  typedef Dune::FieldVector< double , 1 > ValueType;

  {
    Dune::Fem::DynamicArray< Dune::Fem::DynamicArray< ValueType > > vector;
    vector.setMemoryFactor( 1.1 );

    checkVector( vector, 10, 4 );
    checkVector( vector, 7, 6 );
    checkVector( vector, 20, 9 );
    checkVector( vector, 2, 5 );
  }

  {
    Dune::Fem::DynamicArray< ValueType > vector;
    checkPODVector( vector, 10 );
    checkPODVector( vector, 7  );
    checkPODVector( vector, 20 );
    checkPODVector( vector, 2  );
  }

  {
    Dune::Fem::DynamicArray< double > vector;
    checkPODVector( vector, 10 );
    checkPODVector( vector, 7  );
    checkPODVector( vector, 20 );
    checkPODVector( vector, 2  );
  }
  return 0;
}
catch( const Dune::Exception & exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}
