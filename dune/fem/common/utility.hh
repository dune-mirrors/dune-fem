#ifndef DUNE_FEM_COMMON_UTILITY_HH
#define DUNE_FEM_COMMON_UTILITY_HH

#include <algorithm>
#include <type_traits>


namespace Dune
{

  namespace Std
  {

    // sum
    // ---

    template< class T >
    static constexpr T sum ( T a )
    {
      return a;
    }

    template< class T, class ... U >
    static constexpr T sum ( T a, U ... b )
    {
      return a + sum( b ... );
    }


    // sub
    // ---

    template< class T >
    static constexpr T sub ( T a )
    {
      return a;
    }

    template< class T, class ... U >
    static constexpr T sub ( T a, U ... b )
    {
      return a - sum( b ... );
    }


    // max
    // ---

    template< class T >
    static constexpr T max ( T a )
    {
      return a;
    }

    template< class T, class ... U >
    static constexpr T max ( T a, U ... b )
    {
      return std::max( a, max( b ... ) );
    }


    // are_all_same
    // ------------

    //
    // is true if all types in the parameter pack are the same.
    // similar to std::is_same
    //

    template< class ... T >
    struct are_all_same;

    template< class T >
    struct are_all_same< T > : public std::true_type {};

    template< class U, class V, class ... T >
    struct are_all_same< U, V, T ... >
      : public std::integral_constant< bool, std::is_same< U, V >::value &&are_all_same< V, T ... >::value >
    {};

  } // Std

} //  namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_UTILITY_HH
