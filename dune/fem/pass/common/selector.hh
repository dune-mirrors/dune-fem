#ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
#define DUNE_FEM_PASS_COMMON_SELECTOR_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/common/tupleutility.hh>

namespace
{

  // ElementTuple
  // ------------

  /*
   * \brief A helper class that transforms a number of integers to a tuple.
   *        Result type is:
   *  \code
  Dune::tuple< Dune::integral_constant< int, N1 >, ..., Dune::integral_constant< int, Nk > >
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
            class Seed = Dune::tuple<>
          >
  class ElementTuple
  {
    typedef typename Dune::PushBackTuple< Seed, Dune::integral_constant< int, N1 > >::type AccumulatedType;

  public:
    typedef typename ElementTuple< N2, N3, N4, N5, N6, N7, N8, N9, -1, AccumulatedType >::Type Type;
  };

  template< class Seed >
  struct ElementTuple< -1, -1, -1, -1, -1, -1, -1, -1, -1, Seed >
  {
    typedef Seed Type;
  };

} // namespace



namespace Dune
{

  namespace Fem
  {
    
    // Selector
    // --------

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
    struct Selector
    {
      //! \brief tuple consisting of Dune::integral_constant< int, N_i >
      typedef typename ElementTuple< N1, N2, N3, N4, N5, N6, N7, N8, N9 >::Type Type;

      // \brief number of elements in selector
      static const int size = Dune::tuple_size< Type >::value;

      //! \brief check, whether integer N is contained in selector
      template< int N >
      struct Contains
      {
        static const bool v = ( Dune::ContainsType< Type, Dune::integral_constant< int, N > >::value );
      };

    private:
      // Selector is a mere traits class, forbid construction
      Selector();
    };
 
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
