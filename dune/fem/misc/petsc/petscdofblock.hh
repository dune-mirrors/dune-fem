#ifndef DUNE_FEM_PETSCDOFBLOCK_HH
#define DUNE_FEM_PETSCDOFBLOCK_HH

#include <algorithm>
#include <memory>

#include <dune/fem/storage/envelope.hh>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscvector.hh>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {
    template < class PVector >
    class PetscDofBlock;

    /* ========================================
     * class PetscDofBlock::DofProxy
     * =======================================
     */
    template< class PVector >
    class PetscDofProxy
    {
    public:
      typedef PVector  PetscVectorType;
      typedef typename PetscDofBlock< PetscVectorType >::DofProxy  ThisType;
      typedef typename PetscDofBlock< PetscVectorType >::IndexType IndexType;

      static const unsigned int blockSize = PetscVectorType::blockSize;

      // This is needed to put DofProxies in STL (or STL-like) containers...
      PetscDofProxy () = default;

      // ???
      PetscDofProxy ( PetscScalar s ) {}

      // this is called by a friend
      PetscDofProxy ( PetscVectorType &petscVector, IndexType blockIndex, PetscInt indexInBlock )
      : petscVector_( &petscVector ),
        blockIndex_( blockIndex ),
        indexInBlock_( indexInBlock )
      {}

      // Does not change the value of the DoF which this proxy, but it assigns this proxy instance...
      void assign ( const ThisType &other )
      {
        petscVector_ = other.petscVector_;
        blockIndex_ = other.blockIndex_;
        indexInBlock_ = other.indexInBlock_;
      }

      const ThisType& operator= ( ThisType &other )
      {
        setValue( other.getValue() );
        return *this;
      }

      const ThisType& operator= ( PetscScalar val )
      {
        setValue( val );
        return *this;
      }

      const ThisType& operator*= ( const ThisType& other )
      {
        PetscScalar value = getValue() * other.getValue();
        setValue( value );
        return *this;
      }

      const ThisType& operator*= ( const PetscScalar scalar )
      {
        PetscScalar value = getValue() * scalar ;
        setValue( value );
        return *this;
      }

      const ThisType& operator+= ( const ThisType &other )
      {
        setValue( other.getValue(), ADD_VALUES );
        return *this;
      }

      const ThisType& operator+= ( const PetscScalar scalar )
      {
        setValue( scalar, ADD_VALUES );
        return *this;
      }

      const ThisType& operator-= ( const ThisType &other )
      {
        setValue( -other.getValue(), ADD_VALUES );
        return *this;
      }

      const ThisType& operator-= ( const PetscScalar scalar )
      {
        setValue( -scalar, ADD_VALUES );
        return *this;
      }

      // conversion operators
      operator PetscScalar () const
      {
        return valid() ? getValue() : PetscScalar( 0 );
      }

      bool valid() const { return bool( petscVector_ ); }

    protected:
      PetscScalar getValue () const
      {
        assert( valid() );
        PetscInt index = blockSize * petscVector().mappers().ghostIndex( blockIndex_ ) + indexInBlock_;
        PetscScalar ret;
        ::Dune::Petsc::VecGetValues( *petscVector().getGhostedVector(), 1, &index, &ret );
        return ret;
      }

      void setValue ( const PetscScalar &val, InsertMode mode = INSERT_VALUES )
      {
        PetscInt index = blockSize * petscVector().mappers().ghostIndex( blockIndex_ ) + indexInBlock_;
        ::Dune::Petsc::VecSetValue( *( petscVector().getGhostedVector() ), index, val, mode );
        petscVector().hasBeenModified();
      }

      PetscVectorType& petscVector ()
      {
        assert( petscVector_ );
        return *petscVector_;
      }

      const PetscVectorType& petscVector () const
      {
        assert( petscVector_ );
        return *petscVector_;
      }

      // data fields
      PetscVectorType *petscVector_ = nullptr;
      IndexType blockIndex_ = 0;
      PetscInt indexInBlock_ = 0;
    };

    template< class Scalar, class PVector >
    void axpy ( const Scalar &a, const Scalar &x, PetscDofProxy< PVector > proxy )
    {
      proxy += a*x;
    }




    /* ========================================
     * class PetscDofBlock
     * =======================================
     */
    template< class PVector >
    class PetscDofBlock
    {
      typedef PetscDofBlock< PVector > ThisType;

    public:
      typedef PVector PetscVectorType;
      typedef PetscInt IndexType;

      static const unsigned int blockSize = PetscVectorType::blockSize;

      typedef PetscDofProxy< PVector > DofProxy;
      class DofIterator;

      // Needed so that we can put this class in a Dune::Envelope
      typedef std::pair< PetscVectorType&, IndexType > UnaryConstructorParamType;

      // this is the ctor to be used
      PetscDofBlock ( PetscVectorType &petscVector, IndexType blockIndex )
      : petscVector_( petscVector ),
        blockIndex_( blockIndex )
      {}

      // this is the ctor to be used
      PetscDofBlock ( const PetscDofBlock& other )
      : petscVector_( other.petscVector_ ),
        blockIndex_( other.blockIndex_ )
      {}

      // ..or this one, which is semantically equivalent to the above ctor
      explicit PetscDofBlock ( UnaryConstructorParamType arg )
      : petscVector_( arg.first ),
        blockIndex_( arg.second )
      {}

      const PetscDofBlock& operator*= ( const PetscScalar value )
      {
        for( unsigned int i=0; i<blockSize; ++i )
          (*this)[ i ] *= value ;
        return *this;
      }

      PetscDofBlock () = delete;

      IndexType size() const { return blockSize; }

      DofProxy operator [] ( unsigned int index )
      {
        assert( index < blockSize );
        return DofProxy( petscVector_, blockIndex_, index );
      }

      const DofProxy operator [] ( unsigned int index ) const
      {
        assert( index < blockSize );
        return DofProxy( petscVector_, blockIndex_, index );
      }

    private:
      PetscVectorType &petscVector_;
      IndexType blockIndex_;
    };

    //! proper implementation for InStreamInterface of the DofProxy
    //template< class Traits,  class PVector >
    template< class Traits,  class PVector >
    inline OutStreamInterface< Traits >&
      operator<< ( OutStreamInterface< Traits > &out,
                   const PetscDofProxy< PVector >& value )
    {
      // write to stream
      out << PetscScalar( value );
      return out;
    }

    //! proper implementation for InStreamInterface of the DofProxy
    //template< class Traits,  class PVector >
    template< class Traits,  class PVector >
    inline InStreamInterface< Traits >&
      operator>> ( InStreamInterface< Traits > &in,
                   PetscDofProxy< PVector > value )
    {
      PetscScalar val;
      // get value from stream
      in >> val;
      // write value to discrete function
      value = val;
      return in;
    }


    /*
     * This is almost a bidirectional iterator but does not completely satisfy the required
     * interface (see the C++2003 standard, 24.1.4) [no default ctor, no operator->].
     */

    /* ========================================
     * class PetscDofBlock::DofIterator
     * =======================================
     */
    template< class PVector >
    class PetscDofBlock< PVector >::DofIterator
    {
      typedef typename PetscDofBlock< PVector >::DofIterator ThisType;
      typedef PetscDofBlock< PVector > DofBlockType;

      // TODO: get rid of this! we don't like shared pointers. Own a real instance instead!
      typedef std::shared_ptr< DofBlockType > DofBlockSharedPointer;

    public:
      typedef std::input_iterator_tag iterator_category;
      typedef typename DofBlockType::DofProxy value_type;
      typedef std::ptrdiff_t difference_type;
      typedef PetscScalar* pointer;
      typedef PetscScalar& reference;

      typedef PVector PetscVectorType;

      // standard ctor
      DofIterator ( PetscVectorType &petscVector, unsigned int blockIndex, PetscInt indexInBlock = 0 )
      : petscVector_( petscVector ),
        blockIndex_( blockIndex ),
        indexInBlock_( indexInBlock )
      {
        // blockIndex == size denotes the end iterator
        assert( static_cast< std::size_t >( blockIndex ) <= petscVector_.mappers().size() );

        // Is this not the end iterator?
        if( static_cast< std::size_t >( blockIndex ) < petscVector_.mappers().size() )
        {
          resetBlockPtr();
        }
      }

      bool operator== ( const ThisType &other ) const
      {
        return (blockIndex_ == other.blockIndex_) && (indexInBlock_ == other.indexInBlock_);
      }

      bool operator!= ( const ThisType &other ) const { return !this->operator==( other ); };

      value_type operator* () { return block()[ indexInBlock_]; }
      const value_type operator* () const { return block()[ indexInBlock_ ]; }

      // prefix increment
      ThisType& operator++ () { increment(); return *this; }

      // prefix decrement
      ThisType& operator-- () { decrement(); return *this; }

      DofIterator () = delete;

    private:
      PetscVectorType& petscVector () { return petscVector_; }

      void increment ()
      {
        ++indexInBlock_;
        if( static_cast< std::size_t >( indexInBlock_ ) >= DofBlockType::blockSize )
        {
          ++blockIndex_;
          indexInBlock_ = 0;
          if( static_cast< std::size_t >( blockIndex_ ) < petscVector().mappers().size() )
            resetBlockPtr();
        }
      }

      void decrement ()
      {
        assert( blockIndex_ > 0 );
        if( indexInBlock_ == 0 )
        {
          indexInBlock_ = DofBlockType::blockSize - 1;
          --blockIndex_;
          resetBlockPtr();
        }
        else
        {
          --blockIndex_;
        }
      }

      // TODO: iterator should own the block - _no_ new...
      void resetBlockPtr () { blockPtr_.reset( new DofBlockType( *petscVector().block( blockIndex_ ) ) ); }

      DofBlockType& block () const { return *blockPtr_.get(); }

      PetscVectorType &petscVector_;
      unsigned int blockIndex_;
      PetscInt indexInBlock_;
      DofBlockSharedPointer blockPtr_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCDOFBLOCK_HH
