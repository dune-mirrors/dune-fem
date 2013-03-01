#ifndef DUNE_FEM_PASS_COMMON_POINTERTUPLE_HH
#define DUNE_FEM_PASS_COMMON_POINTERTUPLE_HH

#include <dune/common/static_assert.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

#include "tupletypetraits.hh"
#include "tupleutility.hh"

namespace
{

  // DereferenceTuple
  // ----------------

  template< class Tuple,
            class Seed = Dune::tuple<>,
            int index = 0,
            int size = Dune::tuple_size< Tuple >::value
          >
  class DereferenceTuple
  {
    typedef typename Dune::TypeTraits< typename Dune::tuple_element< index, Tuple >::type >::PointeeType & AppendType;
    typedef typename Dune::PushBackTuple< Seed, AppendType >::type AccumulatedType;
    typedef DereferenceTuple< Tuple, AccumulatedType, (index+1), size > NextType;

  public:
    typedef typename Dune::ReferenceTuple< Tuple >::Type Type;

    static Type apply ( Tuple &tuple )
    {
      Seed seed;
      return append( tuple, seed );
    }

  protected:
    template< class, class, int, int > friend class DereferenceTuple;

    static Type append ( Tuple &tuple, Seed &seed )
    {
      typename Dune::tuple_element< index, Tuple >::type pointer = Dune::get< index >( tuple );
      AppendType append = *pointer;
      AccumulatedType next = Dune::tuple_push_back< AppendType >( seed, append );
      return NextType::append( tuple, next );
    }
  };

  template< class Tuple,
            class Seed,
            int size
          >
  struct DereferenceTuple< Tuple, Seed, size, size >
  {
    typedef typename Dune::ReferenceTuple< Tuple >::Type Type;

    static Type apply ( Tuple & )
    {
      return Type();
    }

  protected:
    template< class, class, int, int > friend class DereferenceTuple;

    static Type append ( Tuple &tuple, Seed &seed )
    {
      return seed;
    }
  };

} // namespace



namespace Dune
{

  namespace Fem
  {

    // PointerTuple
    // ------------

    /*
     * \brief Please doc me.
     */
    template< class Tuple >
    class PointerTuple
    {
      typedef PointerTuple< Tuple > ThisType;

      dune_static_assert( TupleTypeTraits< Tuple >::isPointerTuple,
                           "Can not wrap non-pointer tuple." );

   public:
      //! \brief type of pointer tuple
      typedef Tuple PointerType;

      //! \brief type of element
      typedef typename TupleTypeTraits< PointerType >::PointeeTupleType ElementType;

      //! \brief please doc me
      typedef typename DereferenceTuple< PointerType >::Type ReferenceType;

    private:
      // check whether a tuple pointer can be dereferenced
      struct SaneTuplePointer
      {
        SaneTuplePointer () : v_( true ) {}

        template< class Ptr >
        void visit ( Ptr &ptr )
        {
          v_ &= bool( ptr );
        }

        operator bool() const { return v_; }

      private:
        bool v_;
      };

    public:
      //! \brief store null pointer tuple
      PointerTuple ()
      : tuple_( nullptrTuple() )
      {}

      //! \brief store tuple
      explicit PointerTuple ( const Tuple &tuple )
      : tuple_( tuple )
      {}

      //! \brief assignment operator from new tuple
      PointerTuple &operator= ( const ThisType &other )
      {
        tuple_ = other.tuple_;
        return *this;
      }

      //! \brief assignment operator from new tuple
      PointerTuple &operator= ( const Tuple &tuple )
      {
        tuple_ = tuple;
        return *this;
      }

      //! \brief create null pointer
      ReferenceType operator* ()
      {
        return DereferenceTuple< PointerType >::apply( tuple_ );;
      }

      //! \brief create null pointer
      const ReferenceType operator* () const
      {
        return DereferenceTuple< PointerType >::apply( tuple_ );
      }

      //! \brief return true, if internal tuple can be dereferenced
      operator bool () const
      {
        Dune::ForEachValue< PointerType > forEach( tuple_ );
        SaneTuplePointer check;
        forEach.apply( check );
        return check;
      }

    protected:
      static PointerType nullptrTuple ()
      {
        return Dune::NullPointerInitialiser< PointerType >::apply();
      }

    private:
      PointerType tuple_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_POINTERTUPLE_HH
