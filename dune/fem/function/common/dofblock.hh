#ifndef DUNE_FEM_DOFBLOCK_HH
#define DUNE_FEM_DOFBLOCK_HH

#include <utility>
#include <type_traits>

#include <dune/fem/storage/envelope.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class   DofBlockProxy
     *  \ingroup DiscreteFunction
     *  \brief   DoF block proxy for discrete functions with random access to DoFs
     *
     *  Discrete functions are required to support random access to DoF blocks.
     *  If random access to the DoFs is available, DofBlockProxy can be used to
     *  fake the required DoF blocks.
     *
     *  A discrete function willing to use this proxy has to implement a method
     *  called "dof" (and it's constant counterpart) that provides access to the
     *  DoFs via an index.
     *
     *  Use of the DofBlockProxy is shown in VectorDiscreteFunction.
     *
     *  \param  DiscreteFunction  type of the discrete function (for constant
     *                            DoF blocks, use the const version here)
     *  \param  Dof               type of the degrees of freedom (for constant
     *                            DoF blocks, use the const version here)
     *  \param  Size              size of the DoF blocks (block size returned
     *                            from block mapper)
     */
    template< class DiscreteFunction, class Dof, unsigned int Size >
    class DofBlockProxy
    {
      typedef DofBlockProxy< DiscreteFunction, Dof, Size > ThisType;

      friend class Fem :: Envelope< ThisType >;

      typedef DofBlockProxy< const DiscreteFunction, const Dof, Size > ConstDofBlockProxy;
      typedef DofBlockProxy< typename std::remove_const< DiscreteFunction >::type, typename std::remove_const< Dof >::type, Size > NonConstDofBlockProxy;

    public:
      typedef DiscreteFunction DiscreteFunctionType;

      typedef Dof DofType;

      static const unsigned int size = Size;

      typedef unsigned int size_type;

      typedef std::pair< DiscreteFunctionType *, size_type > KeyType;

    protected:
      DofBlockProxy ( const KeyType &key )
      : discreteFunction_( *(key.first) ),
        first_( size * key.second )
      {}

      DofBlockProxy ( const ConstDofBlockProxy &other )
      : discreteFunction_( other.discreteFunction_ ),
        first_( other.first_ )
      {}

      DofBlockProxy ( const NonConstDofBlockProxy &other )
      : discreteFunction_( other.discreteFunction_ ),
        first_( other.first_ )
      {}

    public:
      const ThisType &operator= ( const ConstDofBlockProxy &other )
      {
        for( size_type i = 0; i < size; ++i )
          (*this)[ i ] = other[ i ];
        return *this;
      }

      const ThisType &operator= ( const NonConstDofBlockProxy &other )
      {
        for( size_type i = 0; i < size; ++i )
          (*this)[ i ] = other[ i ];
        return *this;
      }

      const DofType &operator[] ( size_type index ) const
      {
        return discreteFunction_.dof( first_ + index );
      }

      DofType &operator[] ( size_type index )
      {
        return discreteFunction_.dof( first_ + index );
      }

      size_type dim () const
      {
        return size;
      }

    protected:
      DiscreteFunctionType &discreteFunction_;
      const size_type first_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFBLOCK_HH
