#ifndef DUNE_FEM_COMMON_UTILITY_HH
#define DUNE_FEM_COMMON_UTILITY_HH

#include <algorithm>
#include <type_traits>


namespace Dune
{

  namespace Std
  {

    //
    // Set of operations which can performed for an arbitrary number of arguments.
    // Examples:
    //
    // sum( 5, 6, 12, .... )
    // And( true, true, false, ... )
    //
    // or for constant expressions if i... is an integer sequence:
    //
    // sum( std::tuple_element< i, Tuple >::type::value ... )
    //


    // Arithmetical operations

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
      return a - sub( b ... );
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
      return a > max( b ... )? a : max( b ... );
    }


    // Logical operations

    // Or
    // --

    static constexpr bool Or ( bool a )
    {
      return a;
    }

    template < class ... U >
    static constexpr bool Or ( bool a, U ... b )
    {
      return a || Or( b ... );
    }


    // And
    // ---

    static constexpr bool And ( bool a )
    {
      return a;
    }

    template< class B, class ... U >
    static constexpr bool And ( B a, U ... b )
    {
      return a && And( b... );
    }



    // are_all_same
    // ------------

    //
    // is true_type if all types in the parameter pack are the same.
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
