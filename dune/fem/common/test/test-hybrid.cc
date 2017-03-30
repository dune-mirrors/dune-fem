#include <config.h>

#include <iostream>
#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/common/classname.hh>

#include <dune/common/indices.hh>

#include <dune/fem/common/hybrid.hh>

template< class IndexRange >
typename std::decay_t< IndexRange >::Index sum ( IndexRange &&indexRange )
{
  typename std::decay_t< IndexRange >::Index sum( 0 );
  Dune::Hybrid::forEach( std::forward< IndexRange >( indexRange ), [ &sum ] ( auto i ) { sum += i; } );
  return sum;
}

int main ( int argc, char *argv[] )
{
  using Dune::Hybrid::IndexRange;
  using Dune::Hybrid::CompositeIndexRange;

  typedef IndexRange< int, 15 > Range1;
  typedef CompositeIndexRange< CompositeIndexRange< IndexRange< int, 1 >, IndexRange< int, 2 > >, CompositeIndexRange< IndexRange< int, 3 >, IndexRange< int, 4 >, std::make_integer_sequence< int, 5 > > > Range2;

  const auto sum1 = sum( Range1() );
  const auto sum2 = sum( Range2() );
  if( sum1 != sum2 )
  {
    std::cout << "Error: " << sum1 << " != " << sum2 << std::endl;
    return 1;
  }

  typedef std::tuple< std::array< int, 1 >, std::array< int, 2 > > Tuple1;
  typedef std::tuple< std::array< int, 3 >, std::array< int, 4 >, std::array< int, 5 > > Tuple2;
  std::tuple< Tuple1, Tuple2 > tuple;

  using namespace Dune::Indices;
  Dune::Hybrid::forEach( Range2(), [ &tuple ] ( auto i ) {
      std::get< i[ _1 ] >( std::get< i[ _0 ] >( tuple ) )[ i[ _2 ] ] = static_cast< int >( i );
    } );

  return 0;
}
