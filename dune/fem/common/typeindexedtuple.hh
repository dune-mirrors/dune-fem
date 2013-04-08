#ifndef DUNE_FEM_PASS_COMMON_TYPEINDEXEDTUPLE_HH
#define DUNE_FEM_PASS_COMMON_TYPEINDEXEDTUPLE_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

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

  public:
    template< class T >
    struct Value
    {
      typedef typename Dune::tuple_element< Position< T >::value, Tuple >::type Type;
    };

    template< class T >
    struct DUNE_DEPRECATED Get
    : public Value< T >
    {};

    explicit TypeIndexedTuple ( const Tuple &tuple = Tuple() )
    : tuple_( tuple )
    {}

    //! \brief please doc me
    template< class T >
    typename Value< T >::Type &at ()
    {
      return Dune::get< Position< T >::value >( tuple_ );
    }

    //! \brief please doc me
    template< class T >
    const typename Value< T >::Type &at () const
    {
      return Dune::get< Position< T >::value >( tuple_ );
    }

    template< class T >
    DUNE_DEPRECATED typename Value< T >::Type &get () { return at< T >(); }

    template< class T >
    DUNE_DEPRECATED const typename Value< T >::Type &get () const { return at< T >(); }

    //! \brief please doc me
    template< class T >
    typename Value< T >::Type &operator[] ( const T & )
    {
      return at< T >();
    }

    //! \brief please doc me
    template< class T >
    const typename Value< T >::Type &operator[] ( const T & ) const
    {
      return at< T >();
    }

    //! \brief cast to Tuple
    operator Tuple & () { return tuple_; }

    //! \brief cast to const Tuple
    operator const Tuple & () const { return tuple_; }

  private:
    Tuple tuple_;
  };


  // Dune::get for TypeIndexedTuple
  // ------------------------------

  template< int i, class Tuple, class Types >
  typename Dune::tuple_element< i, Tuple >::type &
  get ( TypeIndexedTuple< Tuple, Types > &tuple )
  {
    return Dune::get< i >( static_cast< Tuple & >( tuple ) );
  }

  template< int i, class Tuple, class Types >
  const typename Dune::tuple_element< i, Tuple >::type &
  get ( const TypeIndexedTuple< Tuple, Types > &tuple )
  {
    return Dune::get< i >( static_cast< const Tuple & >( tuple ) );
  }

} // namespace Dune

namespace std { 

  // tuple_element specialization for TypeIndexedTuple
  // -------------------------------------------------

  template< size_t i, class Tuple, class Types >
  struct tuple_element< i, Dune::TypeIndexedTuple< Tuple, Types > > 
  {
    // use types from Tuple, since the get method is specialized to return these 
    typedef typename tuple_element< i, Tuple > :: type type ;
  };
}

#endif // #ifndef DUNE_FEM_PASS_COMMON_TYPEINDEXEDTUPLE_HH
