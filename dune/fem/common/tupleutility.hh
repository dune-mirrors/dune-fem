#ifndef DUNE_FEM_COMMON_TUPLEUTILITY_HH
#define DUNE_FEM_COMMON_TUPLEUTILITY_HH

#include <tuple>
#include <utility>

#include <dune/common/tupleutility.hh>

namespace
{

  // CutOutTuple
  // -----------

  template< class Tuple,
            int begin,
            int length,
            class StartType = std::tuple<>
          >
  class CutOutTuple
  {
    static_assert( (begin+length <= std::tuple_size< Tuple >::value), "Can not cut out tuple of given length" );
    typedef typename Dune::PushBackTuple< StartType, std::tuple_element< begin, Tuple > >::type NextType;

  public:
    typedef typename CutOutTuple< Tuple, (begin+1), (length-1), NextType >::type type;
  };

  template< class Tuple, int begin, class ResultType >
  struct CutOutTuple< Tuple, begin, 0, ResultType >
  {
    typedef ResultType type;
  };

} // namespace



namespace Dune
{

  // PopFrontTuple
  // -------------

  template< class Tuple, int size = std::tuple_size< Tuple >::value >
  struct PopFrontTuple
  {
    static_assert( (size == std::tuple_size< Tuple >::value),
                    "The \"size\" template parameter of PopFrontTuple "
                    "is an implementation detail and should never be "
                    "set explicitly!" );

    typedef typename CutOutTuple< Tuple, 1, (std::tuple_size< Tuple >::value - 1) >::type type;
  };

  template< class Tuple >
  struct PopFrontTuple< Tuple, 0 >
  {
    typedef Tuple type;
  };



  // PopBackTuple
  // ------------

  template< class Tuple, int size = std::tuple_size< Tuple >::value >
  struct PopBackTuple
  {
    static_assert( (size == std::tuple_size< Tuple >::value),
                    "The \"size\" template parameter of PopBackTuple "
                    "is an implementation detail and should never be "
                    "set explicitly!" );

    typedef typename CutOutTuple< Tuple, 0, (std::tuple_size< Tuple >::value - 1) >::type type;
  };

  template< class Tuple >
  struct PopBackTuple< Tuple, 0 >
  {
    typedef Tuple type;
  };



  // tuple_push_back
  // ---------------

  template< typename T, typename... Args >
  inline std::tuple< Args..., T > tuple_push_back ( const std::tuple< Args... > &tup, T t )
  {
    return std::tuple_cat( tup, std::tuple< T > ( t ) );
  }



  // tuple_push_front
  // ----------------

  template< typename T, typename... Args >
  inline std::tuple< T, Args... > tuple_push_front ( const std::tuple< Args... > &tup, T t )
  {
    return std::tuple_cat( std::tuple< T > ( t ), tup );
  }



  // tuple_pop_back
  // --------------

  template< typename Tup, std::size_t... I >
  inline auto tuple_pop_back_impl ( const Tup &tup, const std::index_sequence< I... >& ) ->
    decltype ( std::make_tuple( std::get< I > ( tup )... ) )
  {
    return std::make_tuple ( std::get< I > ( tup )... );
  }

  template< typename T, typename... Args >
  inline auto tuple_pop_back ( const std::tuple< T, Args... > &tup ) ->
    decltype ( tuple_pop_back_impl ( tup , std::make_index_sequence< sizeof...( Args ) > ( ) ) )
  {
    return tuple_pop_back_impl ( tup , std::make_index_sequence< sizeof... ( Args ) > ( ) );
  }



  // tuple_pop_front
  // ---------------

  template< typename Tup, std::size_t... I >
  inline auto tuple_pop_front_impl ( const Tup &tup, const std::index_sequence< I... >& ) ->
    decltype ( std::make_tuple( std::get< I > ( tup )... ) )
  {
    return std::make_tuple ( std::get< I + 1 > ( tup )... );
  }

  template< typename T, typename... Args >
  inline auto tuple_pop_front ( const std::tuple< T, Args... > &tup ) ->
    decltype ( tuple_pop_front_impl ( tup , std::make_index_sequence< sizeof...( Args ) > ( ) ) )
  {
    return tuple_pop_front_impl ( tup , std::make_index_sequence< sizeof... ( Args )  > ( ) );
  }



  // ContainsType
  // ------------

  /*
   * \brief Check wheter given type is contained in tuple.
   *
   * \tparam  Tuple  tuple
   * \tparam  Type   type to search for
   */
  template< class Tuple,
            class Type,
            int N = std::tuple_size< Tuple >::value
          >
  struct ContainsType
  {
    static const bool value = ( std::is_same< typename std::tuple_element< N-1, Tuple >::type, Type >::value
                                || ContainsType< Tuple, Type, N-1 >::value );
  };

  template< class Tuple,
            class Type
          >
  struct ContainsType< Tuple, Type, 0 >
  {
    static const bool value = false;
  };



  // FirstTypeIndexTuple
  // -------------------

  /*
   * \brief Please doc me.
   */
  template< class Tuple,
            class SubTuple,
            class Seed = std::tuple<>,
            int index = 0,
            int size = std::tuple_size< SubTuple >::value
          >
  class FirstTypeIndexTuple
  {
    static_assert( (index == std::tuple_size< Seed >::value),
                    "The \"index\" template parameter of FirstTypeIndexTuple"
                    "is an implementation detail and should never be "
                    "set explicitly!" );

    // get element from selector
    typedef typename std::tuple_element< index, SubTuple >::type Element;
    // find element in pass id tuple
    typedef typename Dune::FirstTypeIndex< Tuple, Element >::type Position;
    // add value to seed
    typedef typename Dune::PushBackTuple< Seed, Position >::type NextSeed;

  public:
    // result type is a tuple of integral constants
    typedef typename FirstTypeIndexTuple< Tuple, SubTuple, NextSeed, (index+1) >::type type;
  };

  template< class Tuple,
            class SubTuple,
            class Seed,
            int size
          >
  struct FirstTypeIndexTuple< Tuple, SubTuple, Seed, size, size >
  {
    typedef Seed type;
  };



  // MakeSubTuple
  // ------------

  /*
   * \brief Please doc me.
   */
  template< class Tuple,
            class Positions,
            class Seed = std::tuple<>,
            int index = 0,
            int size = std::tuple_size< Positions >::value
          >
  class MakeSubTuple
  {
    template< class, class, class, int, int > friend class MakeSubTuple;

    // get pass number for element to append from mapping
    static const int position = std::tuple_element< index, Positions >::type::value;

    // add type to seed
    typedef typename std::tuple_element< position, Tuple >::type AppendType;

    typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;

    typedef MakeSubTuple< Tuple, Positions, AccumulatedType, (index+1), size > NextType;

    static typename NextType::type append ( Tuple &tuple, Seed &seed )
    {
      AppendType append = std::get< position >( tuple );
      AccumulatedType next = tuple_push_back( seed, append );
      return NextType::append( tuple, next );
    }

  public:
    typedef typename NextType::type type;

    static type apply ( Tuple &tuple )
    {
      Seed seed;
      return append( tuple, seed );
    }
  };

  template< class Tuple,
            class Positions,
            class Seed,
            int size >
  class MakeSubTuple< Tuple, Positions, Seed, size, size >
  {
    template< class, class, class, int, int > friend class MakeSubTuple;

    static Seed append ( Tuple &tuple, Seed &seed ) { return seed; }

  public:
    typedef Seed type;

    static type apply ( Tuple & ) { return type(); }
  };



  // TupleToVectorConverter
  // ----------------------

  /** \brief wrapper class to convert a vector of tuples of RangeTypes into something
             that behaves like a vector< RangeType >
  */
  template< class VectorTupleType, int pos >
  class TupleToVectorConverter
  {
  public:
    TupleToVectorConverter ( const TupleToVectorConverter & ) = delete;

    typedef typename VectorTupleType::value_type TupleType;
    typedef typename std::tuple_element< pos, TupleType >::type ValueType;
    typedef ValueType value_type;

    //! constructor
    explicit TupleToVectorConverter ( VectorTupleType &vector )
    : vector_( vector )
    {}

    //! return reference to i-th entry of vector and pos's tuple component
    ValueType &operator [] ( const size_t i )
    {
      using std::get; // for TypeIndexedTuple the function get is in namespace Dune
      assert( i < size() );
      return get< pos >( vector_[ i ] );
    }

    //! return reference to i-th entry of vector and passId's tuple component
    const ValueType &operator [] ( const size_t i ) const
    {
      using std::get; // for TypeIndexedTuple the function get is in namespace Dune
      assert( i < size() );
      return get< pos >( vector_[ i ] );
    }

    //! return size of vector
    size_t size () const { return vector_.size(); }

  protected:
    VectorTupleType &vector_;
  };



  // InstantiateTuple
  // ----------------

  /** \brief Instantiate a tuple of elements with identical, simple
   *         constructors.
   *
   *  \tparam  Tuple  tuple type
   *  \tparam  Key    type of argument to be passed to each constructor
   *  \tparam  Seed   internal template argument
   *  \tparam  len    internal template argument
   *
   * Sample usage:
\code
  Tuple tuple = InstantiateTuple< Tuple, Key >::apply( key );
\endcode
   *  Note that the tuple object is required to be copyable.
   */
  template< class Tuple,
            class Key,
            class Seed = std::tuple<>,
            int len = std::tuple_size< Tuple >::value
          >
  struct InstantiateTuple
  {
    /** \brief create tuple instance
     *
     *  \param[in]  key  argument passed to each constructor
     */
    static Tuple apply ( const Key &key = Key() )
    {
      Seed seed;
      return append( key, seed );
    }

  private:
    template< class, class, class, int > friend struct InstantiateTuple;

    static Tuple append ( const Key &key, Seed &seed )
    {
      static const int index = std::tuple_size< Tuple >::value - len;

      typedef typename std::tuple_element< index, Tuple >::type AppendType;
      typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;

      AccumulatedType next = Dune::tuple_push_back< AppendType >( seed, AppendType( key ) );
      return InstantiateTuple< Tuple, Key, AccumulatedType, len-1 >::append( key, next );
    }
  };

  template< class Tuple, class Key, class Seed >
  struct InstantiateTuple< Tuple, Key, Seed, 0 >
  {
    static Tuple apply ( const Key &key = Key() ) { return Tuple(); }

  private:
    template< class, class, class, int > friend struct InstantiateTuple;

    static Seed append ( const Key &key, Seed &seed ) { return seed; }
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_TUPLEUTILITY_HH
