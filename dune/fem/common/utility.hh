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

    template< class T, std::enable_if_t< std::is_arithmetic< std::decay_t< T > >::value, int > = 0 >
    static constexpr std::decay_t< T > sum ( T a )
    {
      return a;
    }

    template< class T, T a >
    static constexpr std::decay_t< T > sum ( std::integral_constant< T, a > )
    {
      return a;
    }

    template< class T, std::enable_if_t< std::is_enum< std::decay_t< T > >::value, int > = 0 >
    static constexpr std::underlying_type_t< std::decay_t< T > > sum ( T a )
    {
      return a;
    }

    template< class T, class ... U >
    static constexpr auto sum ( T a, U ... b )
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


    // min
    // ---

    template< class T >
    static constexpr T min ( T a )
    {
      return a;
    }

    template< class T, class ... U >
    static constexpr T min ( T a, U ... b )
    {
      return a < min( b ... )? a : min( b ... );
    }


    // Logical operations

    // Or
    // --

    static constexpr bool Or ()
    {
      return false;
    }

    template < class ... U >
    static constexpr bool Or ( bool a, U ... b )
    {
      return a || Or( b ... );
    }


    // And
    // ---

    static constexpr bool And ()
    {
      return true;
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

    // std::is_pod got deprecated

    template<class T>
    struct is_pod
      : public std::integral_constant<bool, std::is_standard_layout<T>::value && std::is_trivial<T>::value>
    {};

  } // Std

  namespace Fem
  {

    namespace detail {

      //! selects Obj::pointSetId if available, otherwise defaultValue (default is -1)
      template <class Obj, int defaultValue = -1 >
      struct SelectPointSetId
      {
      private:
        template <typename T, typename = int>
        struct CheckPointSetId : public std::false_type { };

        template <typename T>
        struct CheckPointSetId<T, decltype((void) T::pointSetId, 0)> : public std::true_type { };

        template <class T, bool>
        struct SelectValue { static const int value = defaultValue; };

        template <class T>
        struct SelectValue< T, true > { static const int value = T::pointSetId;  };
      public:
        static constexpr int value = SelectValue< Obj, CheckPointSetId< Obj >::value >::value;
      };

      //! selects SFS::codegenShapeFunctionSet if available, otherwise defaultValue (default is false)
      template <class SFS, bool defaultValue = false >
      struct IsCodegenShapeFunctionSet
      {
      private:
        template <typename T, typename = int>
        struct CheckCodegenSFS : public std::false_type { };

        template <typename T>
        struct CheckCodegenSFS<T, decltype((void) T::codegenShapeFunctionSet, 0)> : public std::true_type { };

        template <class T, bool>
        struct SelectValue { static const bool value = defaultValue; };

        template <class T>
        struct SelectValue< T, true > { static const bool value = T::codegenShapeFunctionSet;  };
      public:
        static constexpr int value = SelectValue< SFS, CheckCodegenSFS< SFS >::value >::value;
      };
    } // end namespace detail
  } // end namespace Fem

} //  namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_UTILITY_HH
