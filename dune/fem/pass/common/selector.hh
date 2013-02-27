#ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
#define DUNE_FEM_PASS_COMMON_SELECTOR_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

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



  // FindElement
  // -----------

  /*
   * \brief Find element in tuple without throwing an exception
   *        on failure.
   */
  template< class Tuple,
            class Element, 
            int onFailure = -1,
            int start = 0,
            int size = Dune::tuple_size< Tuple >::value
          >
  class FindElement
  {
    static const bool success = ( Dune::is_same< typename Dune::tuple_element< start, Tuple >::type, Element >::value );

  public:
    static const int position = success ? start : FindElement< Tuple, Element, onFailure, start+1, size >::position;
  };

  template< class Tuple,
            class Element,
            int onFailure,
            int size
          >
  struct FindElement< Tuple, Element, onFailure, size, size >
  {
    static const int position = onFailure;
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
        static const bool v = ( FindElement< Type, Dune::integral_constant< int, N > >::position != -1 );
      };

    private:
      // Selector is a mere traits class, forbid construction
      Selector();
    };
 
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_SELECTOR_HH
