#ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
#define DUNE_FEM_PASS_COMMON_SELECTOR_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/selector.hh> instead!"
#endif

#include <type_traits>
#include <tuple>

#include <dune/common/tupleutility.hh>

#include <dune/fem/common/tupleutility.hh>

namespace Dune
{

  namespace Fem
  {

    template< class Seed, int... I >
    struct VariadicSelectorBase;

    /**
     * \brief A helper class that creates a selector tuple from
     *        given pass ids.
     *
     *  \code
     *  std::tuple< std::integral_constant< int, N1 >, ..., std::integral_constant< int, Nk > >
     *  \endcode
     *
     * \note This version does not erase -1 entries.
     */
    template< class... Args, int... I >
    struct VariadicSelectorBase< std::tuple< Args... >, I... >
    {
      //! \brief tuple consisting of std::integral_constant< int, N_i >
      typedef typename JoinTuples< std::tuple<Args...>, std::tuple< std::integral_constant<int,I>...> >::type Type;

      // \brief number of elements in selector
      static const int size = std::tuple_size< Type >::value;

      //! \brief check, whether integer N is contained in selector
      template< int N >
      struct Contains
      {
        static const bool v = ( Dune::ContainsType< Type, std::integral_constant< int, N > >::value );
      };

      VariadicSelectorBase() = delete;
    };

    /**
     * \brief A helper class that creates a selector tuple from
     *        given pass ids.
     */
    template< int... I >
    struct VariadicSelector : public VariadicSelectorBase< std::tuple<>, I... >
    {};



    // ElementTuple
    // ------------
    /*
     * \brief A helper class that transforms a number of integers to a tuple.
     *        Result type is:
     *  \code
     *  std::tuple< std::integral_constant< int, N1 >, ..., std::integral_constant< int, Nk > >
     *  \endcode
     *
     *  \note Terminatory integer values (-1) are cut off from tuple
     */
    template< int N1,
              int N2,
              int N3,
              int N4,
              int N5,
              int N6,
              int N7,
              int N8,
              int N9,
              class Seed = std::tuple<>
            >
    class ElementTuple
    {
      typedef typename Dune::PushBackTuple< Seed, std::integral_constant< int, N1 > >::type AccumulatedType;

    public:
      typedef typename ElementTuple< N2, N3, N4, N5, N6, N7, N8, N9, -1, AccumulatedType >::Type Type;
    };

    template< class Seed >
    class ElementTuple< -1, -1, -1, -1, -1, -1, -1, -1, -1, Seed >
    {
    public:
      typedef Seed Type;
    };

    // Selector
    // --------

    /**
     * \brief A helper class that creates a selector tuple from
     *        given pass ids.
     */
    template< class ElementTupleImp >
    struct SelectorBase;

    template< class ElementTupleImp >
    struct SelectorBase
    {
      //! \brief tuple consisting of std::integral_constant< int, N_i >
      typedef typename ElementTupleImp :: Type Type;

      // \brief number of elements in selector
      static const int size = std::tuple_size< Type >::value;

      //! \brief check, whether integer N is contained in selector
      template< int N >
      struct Contains
      {
        static const bool v = ( Dune::ContainsType< Type, std::integral_constant< int, N > >::value );
      };

      SelectorBase() = delete;
    };

    /**
     * \brief A helper class that creates a selector tuple from
     *        given pass ids.
     */
    template< int N1 = -1,
              int N2 = -1,
              int N3 = -1,
              int N4 = -1,
              int N5 = -1,
              int N6 = -1,
              int N7 = -1,
              int N8 = -1,
              int N9 = -1
            >
    struct Selector : public SelectorBase< Dune::Fem::ElementTuple< N1, N2, N3, N4, N5, N6, N7, N8, N9 > >
    {
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
