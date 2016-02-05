#ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH
#define DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH

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
      void operator() ( const std::size_t index, Value &value ) const
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
      void operator() ( const std::size_t index, Value &value ) const
      {
        axpy( s_, vector_[ index ], value );
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
      void operator() ( const std::size_t index, Value &value ) const
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
      void operator() ( const std::size_t index, Value &value )
      {
        vector_.bind( index, value );
      }

    protected:
      Vector &vector_;
    };

    // DofBlockTraits
    // --------------

    template<class DiscreteFunction >
    struct DofBlockTraits
    {
      typedef typename DiscreteFunction::DofBlockPtrType DofBlockPtrType;
      typedef typename DiscreteFunction::DofBlockType DofBlockType;
    };

    template<class DiscreteFunction >
    struct DofBlockTraits<const DiscreteFunction>
    {
      typedef typename DiscreteFunction::ConstDofBlockPtrType DofBlockPtrType;
      typedef typename DiscreteFunction::ConstDofBlockType DofBlockType;
    };


    // DofBlockFunctor
    // ---------------

    template< class DiscreteFunction, class Functor >
    struct DofBlockFunctor
    {
      typedef typename DofBlockTraits< DiscreteFunction >::DofBlockPtrType DofBlockPtrType;
      typedef typename DofBlockTraits< DiscreteFunction >::DofBlockType DofBlockType;

      static const int blockSize = DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize;

      DofBlockFunctor ( DiscreteFunction &df, Functor &functor )
      : df_( df ), functor_( functor )
      {}

      template < class GlobalKey >
      void operator () ( const std::size_t local, const GlobalKey& globalKey )
      {
        for( int i = 0; i < blockSize; ++i )
          functor_( local*blockSize + i, df_.dofVector()[ globalKey ][ i ] );
      }
    private:
      DiscreteFunction &df_;
      Functor &functor_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_FUNCTOR_HH
