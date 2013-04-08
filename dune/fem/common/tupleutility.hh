#ifndef DUNE_FEM_COMMON_TUPLEUTILITY_HH
#define DUNE_FEM_COMMON_TUPLEUTILITY_HH

#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#ifdef HAVE_TUPLE
#define DUNE_OPEN_TUPLE_NAMESPACE namespace std {
#define DUNE_CLOSE_TUPLE_NAMESPACE }
#elif defined HAVE_TR1_TUPLE // #ifdef HAVE_TUPLE
#define DUNE_OPEN_TUPLE_NAMESPACE namespace std { namespace tr1 {
#define DUNE_CLOSE_TUPLE_NAMESPACE } }
#else // #elif defined HAVE_TR1_TUPLE // #ifdef HAVE_TUPLE
#define DUNE_OPEN_TUPLE_NAMESPACE namespace Dune {
#define DUNE_CLOSE_TUPLE_NAMESPACE }
#endif // #else // #elif defined HAVE_TR1_TUPLE // #ifdef HAVE_TUPLE


namespace
{

  // CutOutTuple
  // -----------

  template< class Tuple,
            int begin,
            int length,
            class StartType = Dune::tuple<>
          >
  class CutOutTuple
  {
    dune_static_assert( (begin+length <= Dune::tuple_size< Tuple >::value),
                        "Can not cut out tuple of given length" );
    typedef typename Dune::PushBackTuple< StartType, Dune::tuple_element< begin, Tuple > >::type NextType;

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

  template< class Tuple, int size = Dune::tuple_size< Tuple >::value >
  struct PopFrontTuple
  {
    dune_static_assert( (size == Dune::tuple_size< Tuple >::value),
                        "The \"size\" template parameter of PopFrontTuple "
                        "is an implementation detail and should never be "
                        "set explicitly!" );

    typedef typename CutOutTuple< Tuple, 1, (Dune::tuple_size< Tuple >::value - 1) >::type type;
  };

  template< class Tuple >
  struct PopFrontTuple< Tuple, 0 >
  {
    typedef Tuple type;
  };



  // PopBackTuple
  // ------------

  template< class Tuple, int size = Dune::tuple_size< Tuple >::value >
  struct PopBackTuple
  {
    dune_static_assert( (size == Dune::tuple_size< Tuple >::value),
                        "The \"size\" template parameter of PopBackTuple "
                        "is an implementation detail and should never be "
                        "set explicitly!" );

    typedef typename CutOutTuple< Tuple, 0, (Dune::tuple_size< Tuple >::value - 1) >::type type;
  };

  template< class Tuple >
  struct PopBackTuple< Tuple, 0 >
  {
    typedef Tuple type;
  };



  // tuple_push_back
  // ---------------

  template< class T1 >
  inline tuple< T1 > tuple_push_back ( const tuple<> &t, T1 t1 )
  {
    return tuple< T1 >( t1 );
  }

  template< class T2, class T1 >
  inline tuple< T1, T2 > tuple_push_back ( const tuple< T1 > &t, T2 t2 )
  {
    return tuple< T1, T2 >( get< 0 >( t ), t2 );
  }

  template< class T3, class T1, class T2 >
  inline tuple< T1, T2, T3 > tuple_push_back ( const tuple< T1, T2 > &t, T3 t3 )
  {
    return tuple< T1, T2, T3 >( get< 0 >( t ), get< 1 >( t ), t3 );
  }

  template< class T4, class T1, class T2, class T3 >
  inline tuple< T1, T2, T3, T4 > tuple_push_back ( const tuple< T1, T2, T3 > &t, T4 t4 )
  {
    return tuple< T1, T2, T3, T4 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), t4 );
  }

  template< class T5, class T1, class T2, class T3, class T4 >
  inline tuple< T1, T2, T3, T4, T5 > tuple_push_back ( const tuple< T1, T2, T3, T4 > &t, T5 t5 )
  {
    return tuple< T1, T2, T3, T4, T5 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), t5 );
  }

  template< class T6, class T1, class T2, class T3, class T4, class T5 >
  inline tuple< T1, T2, T3, T4, T5, T6 > tuple_push_back ( const tuple< T1, T2, T3, T4, T5 > &t, T6 t6 )
  {
    return tuple< T1, T2, T3, T4, T5, T6 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), t6 );
  }

  template< class T7, class T1, class T2, class T3, class T4, class T5, class T6 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7 > tuple_push_back ( const tuple< T1, T2, T3, T4, T5, T6 > &t, T7 t7 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), t7 );
  }

  template< class T8, class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7, T8 > tuple_push_back ( const tuple< T1, T2, T3, T4, T5, T6, T7 > &t, T8 t8 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7, T8 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), t8 );
  }

  template< class T9, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > tuple_push_back ( const tuple< T1, T2, T3, T4, T5, T6, T7, T8 > &t, T9 t9 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), get< 7 >( t ), t9 );
  }



  // tuple_push_front
  // ----------------

  template< class T1 >
  inline tuple< T1 > tuple_push_front ( const tuple<> &t, T1 t1 )
  {
    return tuple< T1 >( t1 );
  }

  template< class T1, class T2 >
  inline tuple< T1, T2 > tuple_push_front ( const tuple< T2 > &t, T1 t1 )
  {
    return tuple< T1, T2 >( t1, get< 0 >( t ) );
  }

  template< class T1, class T2, class T3 >
  inline tuple< T1, T2, T3 > tuple_push_front ( const tuple< T2, T3 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3 >( t1, get< 0 >( t ), get< 1 >( t ) );
  }

  template< class T1, class T2, class T3, class T4 >
  inline tuple< T1, T2, T3, T4 > tuple_push_front ( const tuple< T2, T3, T4 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5 >
  inline tuple< T1, T2, T3, T4, T5 > tuple_push_front ( const tuple< T2, T3, T4, T5 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4, T5 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6 >
  inline tuple< T1, T2, T3, T4, T5, T6 > tuple_push_front ( const tuple< T2, T3, T4, T5, T6 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4, T5, T6 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7 > tuple_push_front ( const tuple< T2, T3, T4, T5, T6, T7 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7, T8 > tuple_push_front ( const tuple< T2, T3, T4, T5, T6, T7, T8 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7, T8 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > tuple_push_front ( const tuple< T2, T3, T4, T5, T6, T7, T8, T9 > &t, T1 t1 )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 >( t1, get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), get< 7 >( t ) );
  }



  // tuple_pop_back
  // ----------------

  template< class T1 >
  inline tuple<> tuple_pop_back ( const tuple< T1 > &t )
  {
    return tuple<>();
  }

  template< class T1, class T2 >
  inline tuple< T1 > tuple_pop_back ( const tuple< T1, T2 > &t )
  {
    return tuple< T1 >( get< 0 >( t ) );
  }

  template< class T1, class T2, class T3 >
  inline tuple< T1, T2 > tuple_pop_back ( const tuple< T1, T2, T3 > &t )
  {
    return tuple< T1, T2 >( get< 0 >( t ), get< 1 >( t ) );
  }

  template< class T1, class T2, class T3, class T4 >
  inline tuple< T1, T2, T3 > tuple_pop_back ( const tuple< T1, T2, T3, T4 > &t )
  {
    return tuple< T1, T2, T3 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5 >
  inline tuple< T1, T2, T3, T4 > tuple_pop_back ( const tuple< T1, T2, T3, T4, T5 > &t )
  {
    return tuple< T1, T2, T3, T4 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6 >
  inline tuple< T1, T2, T3, T4, T5 > tuple_pop_back ( const tuple< T1, T2, T3, T4, T5, T6 > &t )
  {
    return tuple< T1, T2, T3, T4, T5 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
  inline tuple< T1, T2, T3, T4, T5, T6 > tuple_pop_back ( const tuple< T1, T2, T3, T4, T5, T6, T7 > &t )
  {
    return tuple< T1, T2, T3, T4, T5, T6 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7 > tuple_pop_back ( const tuple< T1, T2, T3, T4, T5, T6, T7, T8 > &t )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9 >
  inline tuple< T1, T2, T3, T4, T5, T6, T7, T8 > tuple_pop_back ( const tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > &t )
  {
    return tuple< T1, T2, T3, T4, T5, T6, T7, T8 >( get< 0 >( t ), get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), get< 7 >( t ) );
  }



  // tuple_pop_front
  // ----------------

  template< class T1 >
  inline tuple<> tuple_pop_front ( const tuple< T1 > &t )
  {
    return tuple<>();
  }

  template< class T1, class T2 >
  inline tuple< T2 > tuple_pop_front ( const tuple< T1, T2 > &t )
  {
    return tuple< T2 >( get< 1 >( t ) );
  }

  template< class T1, class T2, class T3 >
  inline tuple< T2, T3 > tuple_pop_front ( const tuple< T1, T2, T3 > &t )
  {
    return tuple< T2, T3 >( get< 1 >( t ), get< 2 >( t ) );
  }

  template< class T1, class T2, class T3, class T4 >
  inline tuple< T2, T3, T4 > tuple_pop_front ( const tuple< T1, T2, T3, T4 > &t )
  {
    return tuple< T2, T3, T4 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5 >
  inline tuple< T2, T3, T4, T5 > tuple_pop_front ( const tuple< T1, T2, T3, T4, T5 > &t )
  {
    return tuple< T2, T3, T4, T5 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6 >
  inline tuple< T2, T3, T4, T5, T6 > tuple_pop_front ( const tuple< T1, T2, T3, T4, T5, T6 > &t )
  {
    return tuple< T2, T3, T4, T5, T6 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
  inline tuple< T2, T3, T4, T5, T6, T7 > tuple_pop_front ( const tuple< T1, T2, T3, T4, T5, T6, T7 > &t )
  {
    return tuple< T2, T3, T4, T5, T6, T7 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
  inline tuple< T2, T3, T4, T5, T6, T7, T8 > tuple_pop_front ( const tuple< T1, T2, T3, T4, T5, T6, T7, T8 > &t )
  {
    return tuple< T2, T3, T4, T5, T6, T7, T8 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), get< 7 >( t ) );
  }

  template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9 >
  inline tuple< T2, T3, T4, T5, T6, T7, T8, T9 > tuple_pop_front ( const tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > &t )
  {
    return tuple< T2, T3, T4, T5, T6, T7, T8, T9 >( get< 1 >( t ), get< 2 >( t ), get< 3 >( t ), get< 4 >( t ), get< 5 >( t ), get< 6 >( t ), get< 7 >( t ), get< 8 >( t ) );
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
            int N = Dune::tuple_size< Tuple >::value
          >
  struct ContainsType
  {
    static const bool value = ( Dune::is_same< typename Dune::tuple_element< N-1, Tuple >::type, Type >::value
                                || ContainsType< Tuple, Type, N-1 >::position );
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
            class Seed = Dune::tuple<>,
            int index = 0,
            int size = Dune::tuple_size< SubTuple >::value
          >
  class FirstTypeIndexTuple
  {
    dune_static_assert( (index == Dune::tuple_size< Seed >::value),
                        "The \"index\" template parameter of FirstTypeIndexTuple"
                        "is an implementation detail and should never be "
                        "set explicitly!" );

    // get element from selector
    typedef typename Dune::tuple_element< index, SubTuple >::type Element;
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
            class Seed = Dune::tuple<>,
            int index = 0,
            int size = Dune::tuple_size< Positions >::value
          >
  class MakeSubTuple
  {
    template< class, class, class, int, int > friend class MakeSubTuple;

    // get pass number for element to append from mapping
    static const int position = Dune::tuple_element< index, Positions >::type::value;

    // add type to seed
    typedef typename Dune::tuple_element< position, Tuple >::type AppendType;

    typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;

    typedef MakeSubTuple< Tuple, Positions, AccumulatedType, (index+1), size > NextType;

    static typename NextType::type append ( Tuple &tuple, Seed &seed )
    {
      AppendType append = Dune::template get< position >( tuple );
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
    // no copying
    TupleToVectorConverter ( const TupleToVectorConverter & );

  public:
    typedef typename VectorTupleType::value_type TupleType;
    typedef typename tuple_element< pos, TupleType >::type ValueType;
    typedef ValueType value_type;

    //! constructor
    explicit TupleToVectorConverter ( VectorTupleType &vector )
    : vector_( vector )
    {}

    //! return reference to i-th entry of vector and pos's tuple component
    ValueType &operator [] ( const size_t i )
    {
      assert( i < size() );
      return Dune::get< pos >( vector_[ i ] );
    }

    //! return reference to i-th entry of vector and passId's tuple component
    const ValueType &operator [] ( const size_t i ) const
    {
      assert( i < size() );
      return get< pos >( vector_[ i ] );
    }

    //! return size of vector 
    size_t size () const { return vector_.size(); }

  protected:
    VectorTupleType &vector_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_TUPLEUTILITY_HH
