#ifndef DUNE_FEM_PASS_COMMON_WRAPPER_HH
#define DUNE_FEM_PASS_COMMON_WRAPPER_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/typetraits.hh>

namespace Dune
{

  namespace Fem
  {

    // Wrapper
    // -------

    /*
     * \brief Wrapper for Dune::tuple objects. Elements are adressed
     *        by pass ids in the selector.
     *
     */
    template< class Tuple, class Selector >
    class Wrapper
    {
      template< int id >
      struct Position
      {
        static const int v = Dune::FirstTypeIndex< Dune::integral_constant< int, id >, Selector >::value; 
      };

    public:
      template< int id >
      struct Get
      {
        typedef typename Dune::tuple_element< Position< id >::v, Tuple >::type Type;
      };

      Wrapper ( Tuple &tuple ) : tuple_( tuple ) {}

      template< int id >
      inline const typename Get< id >::Type &get () const
      {
        return Dune::get< Position< id >::v >( tuple_ );
      }

      template< int id >
      inline typename Get< id >::Type &get ()
      {
        return Dune::get< Position< id >::v >( tuple_ );
      }

      template< int id >
      inline const typename Get< id >::Type &
      operator[] ( const Dune::integral_constant< int, id > & ) const
      {
        return get< id >();
      }

      template< int id >
      inline typename Get< id >::Type &
      operator[] ( const Dune::integral_constant< int, id > & )
      {
        return get< id >();
      }

    private:
      Tuple &tuple_;
    };

  } // namespace Fem

} // namespace Dune


#endif // #ifndef DUNE_FEM_PASS_COMMON_WRAPPER_HH
