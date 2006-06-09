#ifndef DUNE_SELECTION_HH
#define DUNE_SELECTION_HH

#include "tuples.hh"
#include <dune/common/utility.hh> 
#include <dune/common/misc.hh>

namespace Dune {
  /**
   * @brief Converts a selector definition to a set of pairs.
   *
   * A selector is a tuple where the integer constants are converted using
   * Int2Type and the end marker is mapped from -1 to Nil.
   */
  template <int N1, 
            int N2, 
            int N3, 
            int N4, 
            int N5, 
            int N6, 
            int N7, 
            int N8, 
            int N9>
  struct SelectorToPairs {
    typedef Pair<
      Int2Type<N1>, 
      typename SelectorToPairs<N2, N3, N4, N5, N6, N7, N8, N9, -1>::Type
    > Type;
  };

  /*
   * @brief Specialisation for the closure.
   *
   * For selectors, -1 serves as the termination marker.
   */
  template <int N1>
  struct SelectorToPairs<N1, -1, -1, -1, -1, -1, -1, -1, -1> {
    typedef Pair<Int2Type<N1>, Nil> Type;
  };

  /**
   * @brief A list of compile-time constants to select entries of a tuple.
   */
  template <int N1,
            int N2 = -1, 
            int N3 = -1, 
            int N4 = -1, 
            int N5 = -1, 
            int N6 = -1, 
            int N7 = -1, 
            int N8 = -1, 
            int N9 = -1>
  struct Selector : 
    public SelectorToPairs<N1, N2, N3, N4, N5, N6, N7, N8, N9> 
  {
    typedef typename SelectorToPairs<
      N1, N2, N3, N4, N5, N6, N7, N8, N9>::Type Base;
  private:
    //! A Selector contains pure type information, so there is no need to have
    //! instances
    Selector();
  };

  template <class SelectorType>
  struct MaxIndex {
    enum { value = -1 };
  };

  template <class Head, class Tail>
  struct MaxIndex<Pair<Head, Tail> > {
    enum { value = (static_cast<int>(Head::value) > 
                    static_cast<int>(MaxIndex<Tail>::value) ?
                    static_cast<int>(Head::value) : 
                    static_cast<int>(MaxIndex<Tail>::value)) };
  };

  template <>
  struct MaxIndex<Nil> {
    enum { value = -1 };
  };
} // end namespace Dune

#endif
