#ifndef DUNE_FEM_PASS_COMMON_FILTER_HH
#define DUNE_FEM_PASS_COMMON_FILTER_HH

#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#include "tupleutility.hh"
#include "tupletypetraits.hh"

namespace
{

  // PassTree
  // --------

  /*
   * \brief Get pass tree for given pass. Result type is
   *        a tuple consisting of passes beginning from
   *        StartPass.
   */
  template< class Pass,
            int passNum = Pass::passNum,
            class Seed = Dune::tuple<>
          >
  class PassTree 
  {
    dune_static_assert( (passNum == Pass::passNum),
                        "The \"passNum\" template parameter of PassTree"
                        "is an implementation detail and should never be "
                        "set explicitly!" );

    typedef typename Pass::PreviousPassType PreviousPass;

  public:
    typedef typename PassTree< PreviousPass, PreviousPass::passNum,
                               typename Dune::PushFrontTuple< Seed, Pass >::type 
                             >::Type Type;
  };

  // StartPass has pass number 0
  template< class Pass, class Seed >
  struct PassTree< Pass, 0, Seed > 
  {
    typedef typename Dune::PushFrontTuple< Seed, Pass >::type Type;
  };



  // PasIds
  // ------

  /*
   * \brief Extract pass ids from pass tree.
   */
  template< class PassTree >
  class PassIds
  {
    template< class Accumulated, class Pass >
    class AddPassId
    {
      static const int passNum = Pass::passId;

    public:
      // pass ids are inserted consecutively!
      typedef typename Dune::PushBackTuple< Accumulated, Dune::integral_constant< int, passNum > >::type type;
    };

  public:
    typedef typename Dune::ReduceTuple< AddPassId, PassTree >::type Type;
  };



  // Mapping
  // -------

  /*
   * \brief Please doc me.
   */
  template< class PassIds,
            class Selector,
            class Seed = Dune::tuple<>,
            int index = 0,
            int size = Dune::tuple_size< Selector >::value
          >
  class Mapping
  {
    dune_static_assert( (index == Dune::tuple_size< Seed >::value),
                        "The \"index\" template parameter of Mapping"
                        "is an implementation detail and should never be "
                        "set explicitly!" );

    // get element from selector
    typedef typename Dune::tuple_element< index, Selector >::type Element;
    // find element in pass id tuple
    typedef typename Dune::FirstTypeIndex< PassIds, Element >::type Position;
    // add value to seed
    typedef typename Dune::PushBackTuple< Seed, Position >::type NextSeed;

  public:
    // result type is a tuple of pass numbers (wrapped in integral_constant)
    // in the order as definded by the selector
    typedef typename Mapping< PassIds, Selector, NextSeed, (index+1) >::Type Type;
  };

  template< class PassIds,
            class Selector,
            class Seed,
            int size
          >
  struct Mapping< PassIds, Selector, Seed, size, size >
  {
    typedef Seed Type;
  };



  // ExtractType
  // -----------

  /*
   * \brief Please doc me.
   */
  template< class Argument,
            class Mapping,
            class Seed = Dune::tuple<>,
            int index = 0,
            int size = Dune::tuple_size< Mapping >::value
          >
  class ExtractType
  {
    dune_static_assert( ( Dune::tuple_size< Argument >::value >= Dune::tuple_size< Mapping >::value ),
                        "Argument is not long enough" );

    // get pass number for element to append from mapping
    static const int position = Dune::tuple_element< index, Mapping >::type::value;

    // add type to seed
    typedef typename Dune::tuple_element< position, Argument >::type AppendType;
    typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;

    typedef ExtractType< Argument, Mapping, AccumulatedType, (index+1), size > NextType;

  public:
    typedef typename NextType::Type Type;

    static Type apply ( Argument &argument )
    {
      Seed seed;
      return append( argument, seed );
    }

  protected:
    template< class, class, class, int, int > friend class ExtractType;

    static Type append ( Argument &argument, Seed &seed )
    {
      AppendType append = Dune::template get< position >( argument );
      AccumulatedType next = tuple_push_back( seed, append );
      return NextType::append( argument, next );
    }
  };

  template< class Argument,
            class Mapping,
            class ResultType,
            int size >
  struct ExtractType< Argument, Mapping, ResultType, size, size >
  {
    typedef ResultType Type;

  protected:
    template< class, class, class, int, int > friend class ExtractType;

    static Type append ( Argument &argument, ResultType &result )
    {
      return result;
    }
  };

} // namespace



namespace Dune
{

  namespace Fem
  {

    // Filter
    // ------

    /**
     * \brief Filter for total argument tuple for given pass and selector.
     *        Filter::ResultType is a tuple of discrete functions.
     *
     * \tparam  Argument       total argument tuple (see dune/fem/pass/pass.hh)
     * \tparam  Pass           pass type
     * \tparam  SelectorTuple  selector tuple
     *
     */
    template< class Argument, class Pass, class SelectorTuple >
    struct Filter
    {
      //! \brief type of complete argument list
      typedef Argument ArgumentType;
      //! \brief pass type
      typedef Pass PassType;
      //! \brief tuple of pass ids
      typedef SelectorTuple Selector;

    private:
      // get pass tree
      typedef typename PassTree< PassType >::Type PassTreeType;

    public:
      /** \brief a tuple that describes the mapping of local indicies
       *         (in ResultType) to pass numbers (wrapped in
       *         integral_constant)
       */
      typedef typename ::Mapping< typename PassIds< PassTreeType >::Type, Selector >::Type Mapping;

      //! \brief result type is a tuple of pointers to discrete functions
      typedef typename ExtractType< ArgumentType, Mapping >::Type ResultType;

      //! \brief please doc me
      static ResultType apply ( ArgumentType &argument )
      {
        return ExtractType< ArgumentType, Mapping >::apply( argument );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_FILTER_HH
