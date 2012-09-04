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

      template< class GlobalKey >
      void operator() ( const std::size_t local, const GlobalKey &globalKey )
      {
        assign_( globalKey, array_[ local ] );
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
      : localFixed_( i ),
        value_( value )
      {}

      void operator() ( const std::size_t local, const Value &globalKey )
      {
        if( local == localFixed_ )
          value_ = globalKey;
      }

    private:
      std::size_t localFixed_;
      Value &value_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_FUNCTOR_HH
