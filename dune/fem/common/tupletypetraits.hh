#ifndef DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH
#define DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH

#include <tuple>
#include <type_traits>

#include <dune/common/hybridutilities.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/common/tupleutility.hh>

namespace
{

  // CheckAllElements
  // ----------------

  template< class Tuple,
            template< class > class Predicate,
            int N = std::tuple_size< Tuple >::value
          >
  struct CheckAllElements
  {
    static const bool value = ( Predicate< typename std::tuple_element< N-1, Tuple >::type >::value
                                && CheckAllElements< Tuple, Predicate, (N-1) >::value );
  };

  template< class Tuple,
            template< class > class Predicate
          >
  struct CheckAllElements< Tuple, Predicate, 0 >
  {
    static const bool value = true;
  };

} // namespace



namespace Dune
{

  // SingleTypeTuple
  // ---------------

  /*
   * \brief Make tuple of given length and single type.
   */
  template< class T,
            int size,
            class Seed = std::tuple<>,
            int index = 0
          >
  struct SingleTypeTuple
  {
    typedef typename SingleTypeTuple< T, size, typename Dune::PushBackTuple< Seed, T >::type, (index+1) >::Type Type;
  };

  template< class T,
            int size,
            class Seed
          >
  struct SingleTypeTuple< T, size, Seed, size >
  {
    typedef Seed Type;
  };



  // EmptyTuple
  // ----------

  /*
   * \brief a tuple consisting of Dune::Empty
   */
  template< int size >
  struct EmptyTuple
  : public SingleTypeTuple< Dune::Empty, size >
  {};



  // TupleTypeTraits
  // ---------------

  /*
   * \brief Mimicks Dune::TypeTraits (see dune/common/typetraits.hh) for tuples.
   */
  template< class Tuple >
  class TupleTypeTraits
  {
    template< class T >
    struct IsPointer
    {
      static const bool value = std::is_pointer<T>::value;
    };

    template< class T >
    struct IsReference
    {
      static const bool value = std::is_lvalue_reference<T>::value;
    };

    template< class T >
    struct PointeeTypeEvaluator
    {
      typedef typename std::remove_pointer<T>::type Type;
    };

    template< class T >
    struct ReferredTypeEvaluator
    {
      typedef typename std::remove_reference<T>::type Type;
    };

  public:
    static const bool isPointerTuple = CheckAllElements< Tuple, IsPointer >::value;

    typedef typename std::conditional< isPointerTuple,
                                        typename Dune::ForEachType< PointeeTypeEvaluator, Tuple >::Type,
                                        EmptyTuple< std::tuple_size< Tuple >::value >
                                      >::type PointeeTupleType;

    static const bool isReferenceTuple = CheckAllElements< Tuple, IsReference >::value;

    typedef typename Dune::ForEachType< ReferredTypeEvaluator, Tuple >::Type ReferredTupleType;
  };



  // PointerTuple
  // ------------

  /**
   * \brief Convert a tuple to a tuple of pointer types.
   *
   * \tparam  Tuple  none of the types in Tuple should be of pointer type
   */
  template< class Tuple >
  class PointerTuple
  {
    template< class T >
    struct PointerEvaluator
    {
      typedef typename std::remove_pointer<T>::type * Type;
    };

  public:
    typedef typename Dune::ForEachType< PointerEvaluator, Tuple >::Type Type;
  };



  // ReferenceTuple
  // --------------

  /**
   * \brief Convert a tuple to a tuple of references.
   *
   * \tparam  Tuple  tuple to convert
   */
  template< class Tuple >
  class ReferenceTuple
  {
    template< class T >
    struct ReferenceEvaluator
    {
      typedef typename std::remove_reference<T>::type & Type;
    };

  public:
    typedef typename Dune::ForEachType< ReferenceEvaluator, Tuple >::Type Type;
  };



  // ConstTuple
  // ----------

  /**
   * \brief Add const qualifier to all tuple elements.
   *
   * \tparam  Tuple  tuple to convert
   */
  template< class Tuple >
  class ConstTuple
  {
    template< class T >
    struct ConstEvaluator
    {
      typedef typename std::add_const<typename std::remove_cv<T>::type>::type Type;
    };

  public:
    typedef typename Dune::ForEachType< ConstEvaluator, Tuple >::Type Type;
  };



  // RemoveConstTuple
  // ----------------

  /**
   * \brief Remove const qualifiers from tuple.
   */
  template< class Tuple >
  class RemoveConstTuple
  {
    template< class T >
    struct RemoveConstEvaluator
    {
      typedef typename std::remove_cv< T >::type Type;
    };

  public:
    typedef typename Dune::ForEachType< RemoveConstEvaluator, Tuple >::Type Type;
  };



  // tuple_remove_const
  // ------------------

  /*
   * \brief Convenience structure mimicking Dune::remove_const.
   */
  template< class Tuple >
  struct tuple_remove_const
  {
    typedef typename RemoveConstTuple< Tuple >::Type type;
  };



  // ValidPointerTupleCheck
  // ----------------------

  /**
   * \brief Check whether a pointer tuple can be dereferenced.
   *
   * \tparam  Tuple  tuple of pointer types
   */
  template< class Tuple >
  struct ValidPointerTupleCheck
  {
    static_assert( TupleTypeTraits< Tuple >::isPointerTuple, "Can not check non-pointer tuple." );

    static bool apply ( const Tuple &tuple )
    {
      bool check(true);
      Hybrid::forEach( tuple, [ & ]( auto&& ti ){ check &= static_cast< bool >( ti ); } );
      return check;
    }
  };



  // DereferenceTuple
  // ----------------

  /**
   * \brief Dereference pointer tuple.
   */
  template< class Tuple,
            class Seed = std::tuple<>,
            int index = 0,
            int size = std::tuple_size< Tuple >::value
          >
  class DereferenceTuple
  {
    template< class, class, int, int > friend class DereferenceTuple;

    typedef typename std::remove_pointer< typename std::tuple_element< index, Tuple >::type >::type & AppendType;
    typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;

    typedef DereferenceTuple< Tuple, AccumulatedType, (index+1), size > NextType;

  public:
    typedef typename Dune::ReferenceTuple<
        typename Dune::TupleTypeTraits< Tuple >::PointeeTupleType
      >::Type Type;

    static Type apply ( Tuple &tuple )
    {
      Seed seed;
      return append( tuple, seed );
    }

  private:
    static Type append ( Tuple &tuple, Seed &seed )
    {
      typename std::tuple_element< index, Tuple >::type pointer = std::get< index >( tuple );
      AppendType append = *pointer;
      AccumulatedType next = Dune::tuple_push_back< AppendType >( seed, append );
      return NextType::append( tuple, next );
    }
  };

  template< class Tuple,
            class Seed,
            int size
          >
  class DereferenceTuple< Tuple, Seed, size, size >
  {
    template< class, class, int, int > friend class DereferenceTuple;

  public:
    typedef typename Dune::ReferenceTuple<
        typename Dune::TupleTypeTraits< Tuple >::PointeeTupleType
      >::Type Type;

    static_assert( (std::is_same< Seed, Type >::value), "Failed to dereference pointer tuple." );

    static Type apply ( Tuple & )
    {
      return Type();
    }

  private:
    static Seed append ( Tuple &tuple, Seed &seed ) { return seed; }
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH
