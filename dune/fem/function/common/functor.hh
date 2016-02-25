#ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH
#define DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH

#include <utility>

#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>

namespace Dune
{

  namespace Fem
  {

    // LeftAdd
    // -------

    template< class Vector >
    struct LeftAdd
    {
      LeftAdd ( const Vector &vector )
      : vector_( vector )
      {}

      template< class Value >
      void operator() ( const std::size_t index, Value &&value ) const
      {
        value += vector_[ index ];
      }
    private:
      const Vector &vector_;
    };


    // LeftAddScaled
    // -------------

    template< class Vector, class Scalar >
    struct LeftAddScaled
    {
      LeftAddScaled ( const Vector &vector, const Scalar &s )
      : vector_( vector ),
        s_( s )
      {}

      template< class Value >
      void operator() ( const std::size_t index, Value &&value ) const
      {
        axpy( s_, vector_[ index ], std::forward< Value >( value ) );
      }
    private:
      const Vector &vector_;
      const Scalar &s_;
    };


    // LeftAssign
    // ----------

    template< class Vector >
    struct LeftAssign
    {
      LeftAssign ( const Vector &vector )
      : vector_( vector )
      {}

      template< class Value >
      void operator() ( const std::size_t index, Value &&value ) const
      {
        value = vector_[ index ];
      }
    private:
      const Vector &vector_;
    };


    // AssignReference
    // ---------------

    template< class Vector >
    struct AssignVectorReference
    {
      AssignVectorReference ( Vector &vector )
      : vector_( vector )
      {}

      template< class Value >
      void operator() ( const std::size_t index, Value &&value ) const
      {
        vector_.bind( index, std::forward< Value > ( value ) );
      }

    protected:
      Vector &vector_;
    };


    // DofBlockFunctor
    // ---------------

    template< class DofVector, class Functor >
    struct DofBlockFunctor
    {
      static const int blockSize = DofVector::blockSize;

      DofBlockFunctor ( DofVector &dofVector, Functor functor )
      : dofVector_( dofVector ), functor_( functor )
      {}

      template < class GlobalKey >
      void operator () ( std::size_t local, const GlobalKey& globalKey )
      {
        for( int i = 0; i < blockSize; ++i )
          functor_( local*blockSize + i, dofVector_[ globalKey ][ i ] );
      }
    private:
      DofVector &dofVector_;
      Functor functor_;
    };


    // dofBlockFunctor
    // ---------------

    template< class DofVector, class Functor >
    static inline DofBlockFunctor< DofVector, Functor > dofBlockFunctor ( DofVector &dofVector, Functor functor )
    {
      return DofBlockFunctor< DofVector, Functor >( dofVector, std::move( functor ) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH
