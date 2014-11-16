#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_HH

//- system includes
#include <fstream>
#include <iostream>

//- Dune inlcudes
#include <dune/common/exceptions.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/dofmanager.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#else
#include <dune/fem/space/common/arrays.hh>
#endif

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/dofiterator.hh>
#include <dune/fem/function/localfunction/mutable.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declarations
    template< class DiscreteFunctionSpaceType >
    class ISTLBlockVectorDiscreteFunction;
    template< class DofStorageImp, class DofImp >
    class DofIteratorBlockVectorDiscreteFunction;
    template< class BlockVectorImp, class DofImp >
    class StraightenBlockVector;

    template< class DiscreteFunctionSpaceImp >
    struct DiscreteFunctionTraits< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > >
    {
      enum { localBlockSize = DiscreteFunctionSpaceImp::localBlockSize };
      typedef typename DiscreteFunctionSpaceImp::RangeFieldType RangeFieldType;
      typedef RangeFieldType DofType;

      typedef FieldVector< DofType, localBlockSize > DofBlockType;
      typedef const DofBlockType ConstDofBlockType;
      typedef DofBlockType *DofBlockPtrType;
      typedef const DofBlockType *ConstDofBlockPtrType;

#if HAVE_DUNE_ISTL
      typedef BlockVector< DofBlockType > DofStorageType;
#else
      typedef MutableArray< DofBlockType > DofStorageType;
#endif

      typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;

      //! needs additional mapper because of block structure
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;

      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

      typedef DofIteratorBlockVectorDiscreteFunction< DofStorageType, DofType > DofIteratorType;
      typedef ConstDofIteratorDefault< DofIteratorType > ConstDofIteratorType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;

      typedef StraightenBlockVector< DofStorageType, RangeFieldType > LeakPointerType;
    };

    template< class DofType >
    struct DofTypeWrapper
    {
      enum { blockSize = DofType::dimension };
      static typename DofType::field_type &convert ( DofType &val, int idx )
      {
        return val[ idx ];
      }
    };

    template< >
    struct DofTypeWrapper< double >
    {
      enum { blockSize = 1 };
      static double &convert ( double  &val, int idx ) { return val; }
    };

    template< class BlockVectorImp, class DofImp >
    class StraightenBlockVector
    {
      typedef BlockVectorImp BlockVectorType;
      typedef DofImp DofType;
      typedef typename BlockVectorType::block_type BlockType;
      enum { blockSize = BlockType::dimension };

    public:
      StraightenBlockVector ( BlockVectorImp &vec ) : vec_( vec ) {}

      //! return dof
      DofType &operator[] ( int i )
      {
        assert((i >= 0) && ((unsigned int ) i < blockSize * vec_.size()));
        const unsigned int count = (int) i/blockSize;
        assert( count < vec_.size() );
        const unsigned int local = i%blockSize;
        assert( local < blockSize );
        return vec_[ count ][ local ];
      }

      //! return dof
      const DofType &operator[] ( int i ) const
      {
        assert((i >= 0) && ((unsigned int ) i < blockSize * vec_.size()));
        const unsigned int count = (int) i/blockSize;
        assert( count < vec_.size() );
        const unsigned int local = i%blockSize;
        assert( local < blockSize );
        return vec_[ count ][ local ];
      }

      int size () const
      {
        return vec_.size() * blockSize;
      }

      BlockVectorType &blockVector () { return vec_; }
      const BlockVectorType &blockVector () const { return vec_; }

    protected:
      BlockVectorType &vec_;
    };



    //**********************************************************************
    //! @ingroup BlockVectorDFunction
    //  --ISTLBlockVectorDiscreteFunction
    //
    //! this is one special implementation of a discrete function using an
    //! array for storing the dofs.
    //!
    //**********************************************************************
    template< class DiscreteFunctionSpaceImp >
    class ISTLBlockVectorDiscreteFunction
      : public DiscreteFunctionDefault< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > ThisType;
      typedef DiscreteFunctionDefault< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > >
        BaseType;

      typedef BaseType DiscreteFunctionDefaultType;

    public:
      //! type of discrete functions space
      typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

      //! traits of this type
      typedef DiscreteFunctionTraits< ThisType > Traits;

      //! size of local blocks
      static const int localBlockSize = DiscreteFunctionSpaceType::localBlockSize;

      //! needs additional mapper
      typedef typename Traits::BlockMapperType BlockMapperType;

    private:
      enum { myId_ = 0};

    public:
      //! type of underlying array
      typedef typename Traits::DofStorageType DofStorageType;

      //! my type
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

      //! Type of the range field
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

      /** \brief For ISTL-compatibility */
      typedef typename DofStorageType::block_type block_type;
      typedef typename DofStorageType::block_type::field_type field_type;

      //! Type of the grid
      typedef typename DiscreteFunctionSpaceType::GridType GridType;

      //! the dof iterator type of this function
      typedef typename BaseType::DofIteratorType DofIteratorType;
      typedef typename BaseType::ConstDofIteratorType ConstDofIteratorType;

      //! the type of the unknowns
      typedef typename BaseType::DofType DofType;

      //! type of block stored in block vector
      typedef block_type DofBlockType;

      typedef typename BaseType::ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType::DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType::ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType::LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      //! type of index set
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;

      //! type of LeakPointer
      typedef typename Traits::LeakPointerType LeakPointerType;

    public:
      //! \brief Constructor makes Discrete Function
      ISTLBlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType &f );

      //! \brief Constructor makes Discrete Function with name
      ISTLBlockVectorDiscreteFunction ( const std::string &name,
                                        const DiscreteFunctionSpaceType &f );

      //! \brief Constructor makes Discrete Function
      ISTLBlockVectorDiscreteFunction ( const std::string &name,
                                        const DiscreteFunctionSpaceType &f,
                                        const DofStorageType &data );

      //! \brief Constructor makes Discrete Function from copy
      ISTLBlockVectorDiscreteFunction ( const ThisType &other );

      /**  \brief delete stack of free local functions belonging to this discrete function */
      ~ISTLBlockVectorDiscreteFunction () { delete memObject_; }

      using BaseType::name;
      using BaseType::axpy;
      using BaseType::space;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() */
      DofIteratorType dbegin ();

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() */
      DofIteratorType dend ();

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dbegin() const */
      ConstDofIteratorType dbegin () const;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dend() const */
      ConstDofIteratorType dend () const;

      inline ConstDofBlockPtrType block ( const unsigned int block ) const
      {
        return &(dofVec_[ block ]);
      }

      inline DofBlockPtrType block ( const unsigned int block )
      {
        return &(dofVec_[ block ]);
      }

      inline const RangeFieldType &dof ( const unsigned int index ) const
      {
        return leakPtr_[ index ];
      }

      inline RangeFieldType &dof ( const unsigned int index )
      {
        return leakPtr_[ index ];
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::size */
      inline int size () const
      {
        return dofVec_.size() * localBlockSize;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionDefault::clear */
      void clear();

      /** \copydoc Dune::Fem::DiscreteFunctionDefault::axpy */
      void axpy ( const RangeFieldType &s, const DiscreteFunctionType &g );

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::print */
      void print( std::ostream &out ) const;

      /** \brief return reference to internal block vector
          \return reference to blockVector */
      DofStorageType &blockVector () const { return dofVec_; }

      /** \brief return reference to leak pointer
          \return reference to leakPointer */
      LeakPointerType &leakPointer () { return leakPtr_; }

      /** \brief return const reference to leak pointer
          \return constant reference to leakPointer */
      const LeakPointerType &leakPointer () const { return leakPtr_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression();

    private:
      // allocates dof storage
      DofStorageType &allocateDofStorage();

      typename Traits::LocalDofVectorStackType ldvStack_;

      // single mapper for blocks
      BlockMapperType &blockMapper_;

      // DofStorage that manages the memory for the dofs of this function
      DofStorageInterface *memObject_;

      //! the dofs stored in an array
      DofStorageType &dofVec_;

      //! leak pointer converting block vector to straight vector
      LeakPointerType leakPtr_;
    }; // end class ISTLBlockVectorDiscreteFunction



    //***********************************************************************
    //
    //  --DofIteratorBlockVectorDiscreteFunction
    //
    //***********************************************************************
    /** \brief Iterator over an array of dofs
        \todo Please doc me!
     */
    template< class DofStorageImp, class DofImp >
    class DofIteratorBlockVectorDiscreteFunction
      : public DofIteratorDefault< DofImp, DofIteratorBlockVectorDiscreteFunction< DofStorageImp, DofImp > >
    {
    public:
      typedef DofImp DofType;
      typedef DofStorageImp DofStorageType;
      typedef DofIteratorBlockVectorDiscreteFunction< DofStorageType, DofType > ThisType;

      typedef typename DofStorageType::block_type DofBlockType;
      typedef DofTypeWrapper< DofBlockType > DofWrapperType;
      enum { blockSize = DofWrapperType::blockSize };

      //! Default constructor
      DofIteratorBlockVectorDiscreteFunction ()
        : dofArray_( 0 ),
          count_( 0 ),
          idx_( 0 )
      {}

      //! Constructor (const)
      DofIteratorBlockVectorDiscreteFunction ( const DofStorageType &dofArray, int count )
        : dofArray_( const_cast< DofStorageType * >(&dofArray)),
          count_( count ),
          idx_( 0 ) {}

      //! Constructor
      DofIteratorBlockVectorDiscreteFunction ( DofStorageType &dofArray, int count )
        : dofArray_( &dofArray ),
          count_( count ),
          idx_( 0 ) {}

      //! Copy Constructor
      DofIteratorBlockVectorDiscreteFunction ( const ThisType &other )
        : dofArray_( other.dofArray_ ),
          count_( other.count_ ),
          idx_( other.idx_ )
      {}

      //! Assignment operator
      ThisType &operator= ( const ThisType &other )
      {
        if( &other != this )
        {
          dofArray_ = other.dofArray_;
          count_ = other.count_;
          idx_   = other.idx_;
        }
        return *this;
      }
      //! return dof
      DofType &operator* ()
      {
        assert((count_ >= 0) && (count_ < (unsigned int)dofArray_->size()));
        //return DofWrapperType::convert((*dofArray_)[count_],idx_);
        return ((*dofArray_)[ count_ ][ idx_ ]);
      }

      //! return dof read only
      const DofType &operator* () const
      {
        assert((count_ >= 0) && (count_ < (unsigned int)dofArray_->size()));
        //return DofWrapperType::convert((*dofArray_)[count_],idx_);
        return ((*dofArray_)[ count_ ][ idx_ ]);
      }

      //! go to next dof
      ThisType &operator++ ()
      {
        ++idx_;
        if( idx_ >= blockSize )
        {
          idx_ = 0;
          ++count_;
        }
        return (*this);
      }

      //! compare
      bool operator== ( const ThisType &I ) const
      {
        return (count_ == I.count_) && (idx_ == I.idx_);
      }

      //! compare
      bool operator!= ( const ThisType &I ) const
      {
        return !((*this) == I);
      }

      //! return actual index
      int index () const { return count_; }

      //! set dof iterator back to begin , for const and not const Iterators
      void reset () { count_ = 0; idx_ = 0; }

    private:
      //! the array holding the dofs
      DofStorageType *dofArray_;

      //! index
      mutable size_t count_;
      mutable size_t idx_;

    }; // end DofIteratorBlockVectorDiscreteFunction

    template< class DiscreteFunctionSpaceImp >
    class ManagedDiscreteFunction< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > >
      : public ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > BaseType;

    public:
      typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
      typedef ManagedDiscreteFunction< BaseType > ThisType;
      //! \brief Constructor makes Discrete Function
      ManagedDiscreteFunction ( const DiscreteFunctionSpaceType &f ) : BaseType( f ) {}
      //! \brief Constructor makes Discrete Function with name
      ManagedDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType &f ) : BaseType( name, f ) {}
      //! \brief Constructor makes Discrete Function
      ManagedDiscreteFunction ( const std::string name, const DiscreteFunctionSpaceType &f, const typename BaseType::DofStorageType &data ) : BaseType( name, f, data ) {}
      //! \brief Constructor makes Discrete Function from copy
      ManagedDiscreteFunction ( const ThisType &df ) : BaseType( df ) {}
    };

  } // namespace Fem

} // namespace Dune

#include "blockvectorfunction_inline.hh"
#endif // #ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
