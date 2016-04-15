#ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
#define DUNE_FEM_COMMON_TUPLEFOREACH_HH

#include <tuple>
#include <type_traits>
#include <utility>

namespace Dune
{

  namespace Fem
  {

    template< class Tuple, class Functor, std::size_t ... I >
    void for_each ( Tuple &t, Functor functor, std::index_sequence< I ... > )
    {
      std::ignore = std::make_tuple(
        ( functor( std::get< I >( t ), std::integral_constant< std::size_t, I >() ), I )
        ... );
    }

    template< class Functor, class ... T >
    void for_each ( std::tuple< T ... > &t, Functor functor )
    {
      for_each( t, std::move( functor ), std::index_sequence_for< T ... >() );
    }

  }
}

#endif // #ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
