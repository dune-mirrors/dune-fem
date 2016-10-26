#ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
#define DUNE_FEM_COMMON_TUPLEFOREACH_HH

#include <tuple>

#include "dune/common/hybridutilities.hh"
#include "dune/common/std/utility.hh"

namespace Dune
{
  namespace Fem
  {

    template< typename Functor, typename ... T >
    void for_each( std::tuple< T ... >& t, Functor&& functor )
    {
      Hybrid::forEach( Std::make_index_sequence< sizeof ... ( T ) >{},
        [ & ]( auto i ){ functor( std::get< i >( t ), i ); } );
    }

    template< typename Functor, typename ... T >
    void for_each( const std::tuple< T ... >& t, Functor&& functor )
    {
      Hybrid::forEach( Std::make_index_sequence< sizeof ... ( T ) >{},
        [ & ]( auto i ){ functor( std::get< i >( t ), i ); } );
    }

  }
}

#endif // #ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
