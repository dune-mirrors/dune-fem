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
      template< class T, class U >
      void operator() ( const T &a, U &b ) const
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

      template< class Value >
      void operator() ( const std::size_t local, const Value &value ) const
      {
        assign_( value, array_[ local ] );
      }

    private:
      Array &array_;
      Assign assign_;
    };

    template< class T, class Assign >
    struct AssignFunctor< T *, Assign >
    {
      explicit AssignFunctor ( T *array, const Assign &assign = Assign() )
      : array_( array ),
        assign_( assign )
      {}

      template< class Value >
      void operator() ( const std::size_t local, const Value &value ) const
      {
        assign_( value, array_[ local ] );
      }

    private:
      T *array_;
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

      void operator() ( const std::size_t local, const Value &globalKey ) const
      {
        if( local == localFixed_ )
          value_ = globalKey;
      }

    private:
      std::size_t localFixed_;
      Value &value_;
    };

    template <class Mapper2, class Entity2, class Functor>
    struct MatrixFunctor
    {
      explicit MatrixFunctor( const Mapper2 &mapper2, const Entity2 &entity2, Functor functor )
      : mapper2_(mapper2),
        entity2_(entity2),
        functor_(functor)
      {}
      void operator() ( const std::size_t local1, const typename Functor::GlobalKey &globalKey1 ) const
      {
        functor_.set(local1,globalKey1);
        mapper2_.mapEach(entity2_, functor_);
      }
      private:
      const Mapper2 &mapper2_;
      const Entity2 &entity2_;
      Functor functor_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_FUNCTOR_HH
