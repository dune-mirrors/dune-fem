#ifndef DUNE_FEM_COMMON_TYPEINDEXEDTUPLE_HH
#define DUNE_FEM_COMMON_TYPEINDEXEDTUPLE_HH

#include <tuple>

#include <dune/common/tupleutility.hh>
#include <dune/fem/common/tupleutility.hh>

namespace Dune
{

  // TypeIndexedTuple
  // ----------------

  /*
   * \brief Please doc me.
   */
  template< class Tuple, class Types >
  class TypeIndexedTuple
  {
    template< class T >
    struct Position
    {
      static const int value = Dune::FirstTypeIndex< Types, T >::value;
    };

    // value selector for types that are contained
    template< class T, bool contained >
    struct ValueBase
    {
      typedef typename std::tuple_element< Position< T >::value, Tuple >::type Type;
      static Type& at( Tuple& tuple )
      {
        return std::get< Position< T >::value >( tuple );
      }
      static const Type& at( const Tuple& tuple )
      {
        return std::get< Position< T >::value >( tuple );
      }
    };

    // value selector for types that are not contained
    template< class T >
    struct ValueBase< T, false > // T is not contained in Tuple
    {
      typedef Tuple Type;
      static Type& at( Tuple& tuple )
      {
        return tuple;
      }
      static const Type& at( const Tuple& tuple )
      {
        return tuple;
      }
    };

  public:
    template <class T>
    struct Contains
    {
      static const bool value = ContainsType< Types, T >::value;
    };

    template< class T >
    struct Value : public ValueBase< T, Contains< T >::value >
    {
    };

    explicit TypeIndexedTuple ( const Tuple &tuple = Tuple() )
    : tuple_( tuple )
    {}

    //! \brief return reference to tuple member associated with type T
    template< class T >
    typename Value< T >::Type &at ()
    {
      return Value< T >::at( tuple_ );
    }

    //! \brief return reference to tuple member associated with type T
    template< class T >
    const typename Value< T >::Type &at () const
    {
      return Value< T >::at( tuple_ );
    }

    //! \brief return reference to tuple member associated with type T (integral_constant)
    template< class T >
    typename Value< T >::Type &operator[] ( const T & )
    {
      return at< T >();
    }

    //! \brief return reference to tuple member associated with type T (integral_constant)
    template< class T >
    const typename Value< T >::Type &operator[] ( const T & ) const
    {
      return at< T >();
    }

    //! \brief return true if type T is contained in the tuple
    template< class T >
    bool active( const T& ) const { return Contains< T >::value; }

    //! \brief cast to Tuple
    operator Tuple & () { return tuple_; }

    //! \brief cast to const Tuple
    operator const Tuple & () const { return tuple_; }

  private:
    Tuple tuple_;
  };



  // get for TypeIndexedTuple
  // ------------------------

  template< int i, class Tuple, class Types >
  typename std::tuple_element< i, Tuple >::type &
  get ( Dune::TypeIndexedTuple< Tuple, Types > &tuple )
  {
    return std::get< i >( static_cast< Tuple & >( tuple ) );
  }

  template< int i, class Tuple, class Types >
  const typename std::tuple_element< i, Tuple >::type &
  get ( const Dune::TypeIndexedTuple< Tuple, Types > &tuple )
  {
    return std::get< i >( static_cast< const Tuple & >( tuple ) );
  }

} // namespace Dune


// TODO please check this construction, later.
// At the moment it is needed to make dune-fem-dg compile!
namespace std
{
  template< size_t i, class Tuple, class Types >
  class tuple_element< i, Dune::TypeIndexedTuple< Tuple, Types > >
  {
  public:
    typedef typename std::tuple_element< i, Tuple >::type type;
  };
}

#endif // #ifndef DUNE_FEM_COMMON_TYPEINDEXEDTUPLE_HH
