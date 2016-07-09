#ifndef DUNE_FEM_SPACE_COMMON_DERIVATIVES_HH
#define DUNE_FEM_SPACE_COMMON_DERIVATIVES_HH

#include <type_traits>

namespace Dune
{

  namespace Fem
  {

    namespace Derivatives
    {

      // Type
      // ----

      template< int... order >
      struct Type
      {};



      // contains
      // --------

      template< int order >
      inline static constexpr bool contains ( Type<> ) noexcept
      {
        return false;
      }

      template< int order, int tOrder0, int... tOrder >
      inline static constexpr bool contains ( Type< tOrder0, tOrder... > ) noexcept
      {
        return (order == tOrder0) || contains< order >( Type< tOrder... >() );
      }



      // operator &
      // ----------

      template< int... lOrder >
      inline static constexpr auto operator& ( Type< lOrder... >, Type<> ) noexcept
      {
        return Type< lOrder... >();
      }

      template< int... lOrder, int rOrder0, int... rOrder >
      inline static constexpr auto operator& ( Type< lOrder... >, Type< rOrder0, rOrder... > ) noexcept
      {
        return std::conditional_t< contains< rOrder0 >( Type< lOrder... >() ), Type< lOrder... >, Type< lOrder..., rOrder0 > >() & Type< rOrder... >();
      }



      // operator &
      // ----------

      template< int... lOrder, int... rOrder >
      inline static constexpr auto operator+ ( Type< lOrder... >, Type< rOrder... > ) noexcept
      {
        return Type< lOrder... >() & Type< rOrder... >();
      }



      // aliases
      // -------

      typedef Type< 0 > Values;
      typedef Type< 1 > Jacobians;
      typedef Type< 2 > Hessians;

      namespace
      {

        const Values values = {};
        const Jacobians jacobians = {};
        const Hessians hessians = {};

      } // anonymous namespace

    } // namespace Derivatives

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_DERIVATIVES_HH
