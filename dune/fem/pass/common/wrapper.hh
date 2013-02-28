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
        static const int v = Dune::FirstTypeIndex< Selector, Dune::integral_constant< int, id > >::value; 
      };

    public:
      template< int id >
      struct Get
      {
        static const int position = Position< id >::v;
        typedef typename Dune::tuple_element< position, Tuple >::type Type;
      };

      explicit Wrapper ( const Tuple &tuple ) : tuple_( tuple ) {}

      template< int id >
      inline const typename Get< id >::Type &get () const
      {
        return Dune::get< Position< id >::v >( tuple_ );
      }

      template< int id >
      inline const typename Get< id >::Type &
      operator[] ( const Dune::integral_constant< int, id > & ) const
      {
        return get< id >();
      }

    private:
      const Tuple &tuple_;
    };

  } // namespace Fem

} // namespace Dune


#endif // #ifndef DUNE_FEM_PASS_COMMON_WRAPPER_HH
