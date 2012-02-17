#ifndef DUNE_FEM_MISC_FUNCTOR_HH
#define DUNE_FEM_MISC_FUNCTOR_HH

#include <cstddef>

namespace Dune
{

  namespace Fem
  {

    // DefaultAssign
    // -------------

    struct DefaultAssign
    {
      template< class T >
      void operator() ( const T &a, T &b ) const
      {
        b = a;
      }
    };



    // AssignFunctor
    // -------------

    template< class Array, class Assign = DefaultAssign >
    struct AssignFunctor
    {
      explicit AssignFunctor ( Array &array, const Assign &assign = Assign() )
      : array_( array ),
        assign_( assign )
      {}

      template< class Element >
      void operator() ( const std::size_t i, const Element &e )
      {
        assign( e, array_[ i ] );
      }

    private:
      Array &array_;
      Assign assign_;
    };



    // AssignSingleFunctor
    // -------------------

    template< class Value >
    struct AssignSingleFunctor
    {
      explicit AssignSingleFunctor ( const std::size_t i, Value &value )
      : i_( i ),
        value_( value )
      {}

      void operator() ( const std::size_t i, const Value &v )
      {
        if( i == i_ )
          value_ = v;
      }

    private:
      std::size_t i_;
      Value &value_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_FUNCTOR_HH
