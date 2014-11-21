// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCDOFBLOCK_HH
#define DUNE_FEM_PETSCDOFBLOCK_HH

#include <dune/fem/storage/envelope.hh>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscvector.hh>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {


    // TODO: give this the same interface as Dune::Fem::ReferenceBlockVectorBlock?
    /* ========================================
     * class PetscDofBlock
     * =======================================
     */
    template< typename PVector >
    class PetscDofBlock
    {
      typedef PetscDofBlock< PVector > ThisType;

    public:
      typedef PVector PetscVectorType;
      typedef typename PetscVectorType::PetscDofMappingType PetscDofMappingType;
      typedef typename PetscDofMappingType::IndexType IndexType;

      static const unsigned int blockSize = PetscVectorType::blockSize;

      class DofProxy;
      class DofIterator;

      // Needed so that we can put this class in a Dune::Envelope
      typedef std::pair< PetscVectorType&, IndexType > UnaryConstructorParamType;

      // this is the ctor to be used
      PetscDofBlock ( PetscVectorType &petscVector, IndexType blockIndex )
      : petscVector_( petscVector ),
        blockIndex_( blockIndex )
      {}

      // ..or this one, which is semantically equivalent to the above ctor
      explicit PetscDofBlock ( UnaryConstructorParamType arg )
      : petscVector_( arg.first ),
        blockIndex_( arg.second )
      {}

    private:
      PetscDofBlock ();
      // The standard copy ctor is okay... and needed.

    public:
      DofProxy operator[] ( unsigned int index )
      {
        assert( index < blockSize );
        return DofProxy( petscVector_, blockIndex_, index );
      }

    private:
      /*
       * data fields
       */
      PetscVectorType &petscVector_;
      IndexType blockIndex_;

    };

    /* ========================================
     * class PetscDofBlock::DofProxy
     * =======================================
     */
    template< typename PVector >
    class PetscDofBlock< PVector >::DofProxy
    {
      typedef PVector PetscVectorType;
      typedef typename PetscDofBlock< PetscVectorType >::DofProxy ThisType;
      typedef typename PetscDofBlock< PetscVectorType >::IndexType IndexType;

    public:

      // This is needed to put DofProxies in STL (or STL-like) containers...
      DofProxy ()
      : petscVector_( 0 ),
        blockIndex_( 0 ),
        indexInBlock_( 0 )
      {}

      // this is called by a friend
      DofProxy ( PetscVectorType &petscVector, IndexType blockIndex, PetscInt indexInBlock )
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
      operator PetscScalar () const { return getValue(); }

    private:

      PetscScalar getValue () const
      {
        PetscInt index = blockSize*petscVector().dofMapping().localSlaveMapping( blockIndex_ ) + indexInBlock_;
        PetscScalar ret;
        ::Dune::Petsc::VecGetValues( *petscVector().getGhostedVector(), 1, &index, &ret );
        return ret;
      }

      void setValue ( const PetscScalar &val, InsertMode mode = INSERT_VALUES )
      {
        PetscInt index = blockSize*petscVector().dofMapping().localSlaveMapping( blockIndex_ ) + indexInBlock_;
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
      PetscVectorType *petscVector_;
      IndexType blockIndex_;
      PetscInt indexInBlock_;

    };


    // TODO: to be revised
    //! proper implementation for InStreamInterface of the DofProxy
    //template< class Traits,  class PVector >
    template< class Traits,  class T >
    inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 //const typename PetscDofBlock< PVector >::DofProxy& value )
                 const T& value )
    {
      DUNE_THROW(NotImplemented,"operator>> not implemented for PetscDofBlock< PVector >::DofProxy");
/*
      PetscScalar val;
      in >> val;
      value = val;
*/
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
    template< typename PVector >
    class PetscDofBlock< PVector >::DofIterator
    : public std::iterator< std::input_iterator_tag, PetscScalar >
    {
      typedef typename PetscDofBlock< PVector >::DofIterator ThisType;
      typedef PetscDofBlock< PVector > DofBlockType;

      typedef std::iterator< std::input_iterator_tag, PetscScalar > BaseType;

      // TODO: get rid of this! we don't like shared pointers. Own a real instance instead!
      typedef Dune::shared_ptr< DofBlockType > DofBlockSharedPointer;

    public:
      typedef PVector PetscVectorType;
      typedef typename DofBlockType::DofProxy value_type; // (this overrides the 2nd template parameter of BaseType...)

      // standard ctor
      DofIterator ( PetscVectorType &petscVector, unsigned int blockIndex, PetscInt indexInBlock = 0 )
      : petscVector_( petscVector ),
        blockIndex_( blockIndex ),
        indexInBlock_( indexInBlock )
      {
        // blockIndex == size denotes the end iterator
        assert( blockIndex <= petscVector_.dofMapping().size() );

        // Is this not the end iterator?
        if( blockIndex < petscVector_.dofMapping().size() )
        {
          resetBlockPtr();
        }
      }

      bool operator== ( const ThisType &other ) const
      {
        return blockIndex_ == other.blockIndex_ &&
               indexInBlock_ == other.indexInBlock_;
      }

      bool operator!= ( const ThisType &other ) const { return !this->operator==( other ); };

      value_type operator* () { return block()[ indexInBlock_]; }
      const value_type operator* () const { return block()[ indexInBlock_ ]; }

      // prefix increment
      ThisType& operator++ () { increment(); return *this; }

      // prefix decrement
      ThisType& operator-- () { decrement(); return *this; }

    private:
      // forbidden
      DofIterator ();

      PetscVectorType& petscVector () { return petscVector_; }

      void increment ()
      {
        ++indexInBlock_;
        if( static_cast< unsigned int >( indexInBlock_ ) >= DofBlockType::blockSize )
        {
          ++blockIndex_;
          indexInBlock_ = 0;
          if( static_cast< unsigned int >( blockIndex_ ) < petscVector().dofMapping().size() )
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

      /*
       * data fields
       */
      PetscVectorType &petscVector_;
      unsigned int blockIndex_;
      PetscInt indexInBlock_;
      DofBlockSharedPointer blockPtr_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCDOFBLOCK_HH
