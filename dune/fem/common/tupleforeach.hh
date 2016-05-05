#ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
#define DUNE_FEM_COMMON_TUPLEFOREACH_HH

#include <tuple>

namespace Dune
{
  namespace Fem
  {

    template< typename Functor, std::size_t R >
    struct ForEachImpl
    {
      template< typename ... T >
      ForEachImpl( std::tuple< T ... >& t, Functor&& functor)
      {
        constexpr std::size_t I = sizeof ... ( T ) - R - 1;
        functor( std::get< I >( t ), I );
        ForEachImpl< Functor, R-1 >( t, std::move( functor ) );
      }

      template< typename ... T >
      ForEachImpl( const std::tuple< T ... >& t, Functor&& functor)
      {
        constexpr std::size_t I = sizeof ... ( T ) - R - 1;
        functor( std::get< I >( t ), I );
        ForEachImpl< Functor, R-1 >( t, std::move( functor ) );
      }

    };

    template< typename Functor >
    struct ForEachImpl< Functor, 0 >
    {
      template< typename ... T >
      ForEachImpl( std::tuple< T ... >& t, Functor&& functor )
      {
        constexpr std::size_t I = sizeof ... ( T ) - 1;
        functor( std::get< I > ( t ), I );
      }

      template< typename ... T >
      ForEachImpl( const std::tuple< T ... >& t, Functor&& functor )
      {
        constexpr std::size_t I = sizeof ... ( T ) - 1;
        functor( std::get< I > ( t ), I );
      }
    };

    template< typename Functor, typename ... T >
    void for_each( std::tuple< T ... >& t, Functor functor )
    {
      ForEachImpl< Functor, sizeof ... ( T ) - 1 >( t, std::move( functor ) );
    }

    template< typename Functor, typename ... T >
    void for_each( const std::tuple< T ... >& t, Functor functor )
    {
      ForEachImpl< Functor, sizeof ... ( T ) - 1 >( t, std::move( functor ) );
    }

  }
}

#endif // #ifndef DUNE_FEM_COMMON_TUPLEFOREACH_HH
