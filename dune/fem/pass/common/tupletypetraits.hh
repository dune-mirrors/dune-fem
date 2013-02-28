#ifndef DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH
#define DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH

#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include "tupleutility.hh"

namespace
{

  // CheckAllElements
  // ----------------

  template< class Tuple,
            template< class > class Predicate,
            int N = Dune::tuple_size< Tuple >::value
          >
  struct CheckAllElements
  {
    static const bool value = ( Predicate< typename Dune::tuple_element< N-1, Tuple >::type >::value
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
            class Seed = Dune::tuple<>,
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
      static const bool value = Dune::TypeTraits< T >::isPointer;
    };

    template< class T >
    struct IsReference
    {
      static const bool value = Dune::TypeTraits< T >::isReference;
    };

    template< class T >
    struct PointeeTypeEvaluator
    {
      typedef typename Dune::TypeTraits< T >::PointeeType Type;
    };

    template< class T >
    struct ReferredTypeEvaluator
    {
      typedef typename Dune::TypeTraits< T >::ReferredType Type;
    };

  public:
    static const bool isPointerTuple = CheckAllElements< Tuple, IsPointer >::value;

    typedef typename Dune::SelectType< isPointerTuple,
                                       typename Dune::ForEachType< PointeeTypeEvaluator, Tuple >::Type, 
                                       EmptyTuple< Dune::tuple_size< Tuple >::value > 
                                     >::Type PointeeTupleType;

    static const bool isReferenceTuple = CheckAllElements< Tuple, IsReference >::value;

    typedef typename Dune::ForEachType< ReferredTypeEvaluator, Tuple >::Type ReferredTupleType;
  };



  // tuple_remove_const
  // ------------------

  /*
   * \brief Please doc me.
   */
  template< class Tuple >
  class tuple_remove_const
  {
    template< class T >
    struct RemoveConstEvaluator
    {
      typedef typename Dune::remove_const< T >::type Type;
    };

 public:
    typedef typename Dune::ForEachType< RemoveConstEvaluator, Tuple >::Type type;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_TUPLETYPETRAITS_HH
