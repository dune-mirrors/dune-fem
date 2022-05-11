#ifndef DUNE_FEM_PETSCVECTOR_HH
#define DUNE_FEM_PETSCVECTOR_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/fem/storage/dynamicarray.hh>
#include <dune/fem/storage/envelope.hh>

#include <dune/fem/common/hybrid.hh>

#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofblock.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/petsc.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class DFSpace >
    class PetscVector;


    // External Forward Declarations
    // -----------------------------

    template< class >
    class PetscDofBlock;

    template< class >
    class PetscDofProxy;



    // PetscManagedDofStorage
    // ----------------------

    /*! ManagedDofStorage for PetscDiscreteFunction using PetscVector */
    template < class DiscreteFunctionSpace, class Mapper >
    class PetscManagedDofStorage
      : public ManagedDofStorageImplementation< typename DiscreteFunctionSpace::GridPartType::GridType, Mapper, PetscVector< DiscreteFunctionSpace > >
    {
      typedef typename DiscreteFunctionSpace::GridPartType::GridType GridType;
      typedef Mapper MapperType ;
      typedef PetscVector< DiscreteFunctionSpace > DofArrayType ;
      typedef ManagedDofStorageImplementation< GridType, MapperType, DofArrayType > BaseType;

    public:
      typedef typename BaseType :: SizeType SizeType;

      //! Constructor of ManagedDofStorageImpl, only to call from DofManager
      PetscManagedDofStorage( const DiscreteFunctionSpace& space,
                              const MapperType& mapper )
        : BaseType( space.grid(), mapper, myArray_ ),
          myArray_( space )
      {
      }

      void resize( const bool ) override
      { // do nothing here, only in compress
      }

      void reserve( const SizeType ) override
      {
        // do nothing here, only in compress
      }

      void dofCompress( const bool clearResizedArrays ) override
      {
        myArray_.resize();
        if( clearResizedArrays )
        {
          myArray_.clear();
        }
      }

      //! enable dof compression for dof storage (default is empty)
      void enableDofCompression() override
      {
        std::cerr << "WARNING: PetscVector cannot handle dof compression!" << std::endl;
      }
    protected:
      DofArrayType myArray_;
    };


    // PetscVector
    // -----------

    /*
     * This encapsules a PETSc Vec with ghosts.
     * Some conceptual explanations:
     *  The PETSc vector, as modelled by this class, consists of three parts:
     *  1) the whole PETSc vector, which might be distributed across several MPI processes.
     *     We call this the _global vector_
     *  2) Each process holds a portion of this global vector, we call this part the
     *     _process vector_.
     *  3) And there is a representation of the process vector, which also has 'ghost dof blocks'.
     *     We call this represenation the _ghosted vector_.
     */
    template< class DFSpace >
    class PetscVector
    {
      typedef PetscVector< DFSpace > ThisType;
      friend class PetscDofBlock< ThisType >;
      friend class PetscDofBlock< const ThisType >;
      friend class PetscDofProxy< ThisType >;
      friend class PetscDofProxy< const ThisType >;
    public:
      typedef PetscScalar  value_type ;
      typedef Vec  DofContainerType;

      static constexpr int blockSize = DFSpace::localBlockSize;
      typedef Hybrid::IndexRange< int, blockSize > BlockIndices;

      typedef PetscDofBlock< ThisType >                       DofBlockType;
      typedef PetscDofBlock< const ThisType >                 ConstDofBlockType;
      typedef typename DofBlockType::DofIterator              DofIteratorType;
      typedef typename ConstDofBlockType::DofIterator         ConstDofIteratorType;
      typedef Envelope< DofBlockType >                        DofBlockPtrType;
      typedef Envelope< ConstDofBlockType >                   ConstDofBlockPtrType;
      typedef typename DofBlockType::IndexType                IndexType;

      typedef DofIteratorType       IteratorType;
      typedef ConstDofIteratorType  ConstIteratorType;
      typedef typename DFSpace::RangeFieldType FieldType;
      typedef int  SizeType;

      typedef PetscMappers< DFSpace > MappersType;

      struct FakeDF
      {
        typedef DFSpace DiscreteFunctionSpaceType;
      };

      // need to pass something that contains a typedef DiscreteFunctionSpaceType
      typedef typename DFSpace::template CommDataHandle< FakeDF >::OperationType DefaultCommunicationOperationType;

      // note that Vec is a pointer type so no deep copy is made
      PetscVector ( const DFSpace& space, Vec vec )
        : mappers_( space ), vec_(vec), owner_(false)
      {
        static_assert( DefaultCommunicationOperationType::value == DFCommunicationOperation::copy ||
                       DefaultCommunicationOperationType::value == DFCommunicationOperation::add,
                       "only copy/add are available communication operations for petsc");
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );
      }
      PetscVector ( const DFSpace& space )
        : mappers_( space ), owner_(true)
      {
        static_assert( DefaultCommunicationOperationType::value == DFCommunicationOperation::copy ||
                       DefaultCommunicationOperationType::value == DFCommunicationOperation::add,
                       "only copy/add are available communication operations for petsc");
        // init vector
        init();
      }

      // TODO: think about sequence overflows...
      PetscVector ( const ThisType &other )
        : mappers_( other.mappers_ ), owner_(true)
      {
        // init vector
        init();

        // copy dofs
        assign( other );
      }

      ~PetscVector ()
      {
        if (owner_)
          // destroy vectors
          removeObj();
      }

      std::size_t size () const { return mappers().ghostMapper().size(); }

      void resize( const std::size_t newsize = 0 )
      {
        // TODO: keep old data stored in current vector
        // remove old vectors
        removeObj();

        // initialize new
        init();

        hasBeenModified ();
      }

      void reserve( const std::size_t capacity )
      {
        resize( capacity );
      }

      void hasBeenModified () { ++sequence_; }

      void communicate ()
      {
        communicateFlag_ = true;
      }

      // accessors for the underlying PETSc vectors
      Vec* getVector ()
      {
        communicateIfNecessary();
        return &vec_;
      }

      const Vec* getVector () const
      {
        communicateIfNecessary();
        return &vec_;
      }

      // accessors for the underlying PETSc vectors
      Vec& array ()
      {
        communicateIfNecessary();
        return vec_;
      }

      const Vec& array () const
      {
        communicateIfNecessary();
        return vec_;
      }

      Vec* getGhostedVector ()
      {
        communicateIfNecessary();
        return &ghostedVec_;
      }

      const Vec* getGhostedVector () const
      {
        communicateIfNecessary();
        return &ghostedVec_;
      }

      void beginAssemble()
      {
        ::Dune::Petsc::VecAssemblyBegin( vec_ );
      }

      void endAssemble()
      {
        ::Dune::Petsc::VecAssemblyEnd( vec_ );
      }

      // force communication _now_
      template <class Operation>
      void communicateNow (const Operation& operation) const
      {
        communicateFlag_ = true;
        ++sequence_;
        communicateIfNecessary( operation );
      }

      DofBlockType operator[] ( const IndexType index ) { return DofBlockType( *this,index ); }
      ConstDofBlockType operator[] ( const IndexType index ) const { return ConstDofBlockType( *this,index ); }

      ConstDofBlockPtrType block ( IndexType index ) const { return blockPtr( index ); }
      DofBlockPtrType block ( IndexType index ) { return blockPtr( index ); }

      DofBlockPtrType blockPtr ( IndexType index )
      {
        assert( static_cast< std::size_t >( index ) < mappers().size() );
        return DofBlockPtrType( typename DofBlockType::UnaryConstructorParamType( *this, index ) );
      }

      ConstDofBlockPtrType blockPtr ( IndexType index ) const
      {
        assert( static_cast< std::size_t >( index ) < mappers().size() );
        return ConstDofBlockPtrType( typename ConstDofBlockType::UnaryConstructorParamType( *this, index ) );
      }

      DofIteratorType begin () { return DofIteratorType( *this, 0, 0 ); }
      ConstDofIteratorType begin () const { return ConstDofIteratorType( *this, 0, 0 ); }
      DofIteratorType end () { return DofIteratorType( *this, mappers().size() ); }
      ConstDofIteratorType end () const { return ConstDofIteratorType( *this, mappers().size() ); }

      DofIteratorType beforeEnd ()
      {
        DUNE_THROW(NotImplemented,"PetscVector::beforeEnd is not implemented yet");
        return end();
      }

      ConstDofIteratorType beforeEnd () const
      {
        DUNE_THROW(NotImplemented,"PetscVector::beforeEnd is not implemented yet");
        return end();
      }

      DofIteratorType beforeBegin ()
      {
        DUNE_THROW(NotImplemented,"PetscVector::beforeBegin is not implemented yet");
        return end();
      }
      ConstDofIteratorType beforeBegin () const
      {
        DUNE_THROW(NotImplemented,"PetscVector::beforeBegin is not implemented yet");
        return end();
      }

      DofIteratorType dbegin () { return begin(); }
      ConstDofIteratorType dbegin () const { return begin(); }
      DofIteratorType dend () { return end(); }
      ConstDofIteratorType dend () const { return end(); }

      void clear ()
      {
        ::Dune::Petsc::VecSet( *getVector(), 0. );
        updateGhostRegions();
        vectorIsUpToDateNow();
      }

      PetscScalar operator* ( const ThisType &other ) const
      {
        PetscScalar ret;
        ::Dune::Petsc::VecDot( *getVector(), *other.getVector(), &ret );
        return ret;
      }

      const ThisType& operator+= ( const ThisType &other )
      {
        ::Dune::Petsc::VecAXPY( *getVector(), 1., *other.getVector() );
        updateGhostRegions();
        vectorIsUpToDateNow();
        return *this;
      }

      const ThisType& operator-= ( const ThisType &other )
      {
        ::Dune::Petsc::VecAXPY( *getVector(), -1., *other.getVector() );
        updateGhostRegions();
        vectorIsUpToDateNow();
        return *this;
      }

      const ThisType& operator*= ( PetscScalar scalar )
      {
        ::Dune::Petsc::VecScale( *getVector(), scalar );
        updateGhostRegions();
        vectorIsUpToDateNow();
        return *this;
      }

      const ThisType& operator/= ( PetscScalar scalar )
      {
        assert( scalar != 0 );
        return this->operator*=( 1./scalar );
      }

      void axpy ( const PetscScalar &scalar, const ThisType &other )
      {
        ::Dune::Petsc::VecAXPY( *getVector(), scalar, *other.getVector() );
        hasBeenModified();
      }

      void clearGhost( )
      {
        PetscScalar *array;
        // create local memory
        VecGetArray( ghostedVec_,&array );
        std::fill_n( array + mappers().ghostMapper().interiorSize() * blockSize, mappers().ghostMapper().ghostSize() * blockSize, PetscScalar( 0 ) );
        // destroy memory
        VecRestoreArray( ghostedVec_, &array );
      }

      // debugging; comes in handy to call these 2 methods in gdb
      // doit is only here to prevent the compiler from optimizing these calls away...
      void printGlobal ( bool doit )
      {
          if( !doit )
            return;
          VecView( vec_, PETSC_VIEWER_STDOUT_WORLD );
      }

      void printGhost ( bool doit)
      {
          if( !doit )
            return;

          PetscScalar *array;
          VecGetArray( ghostedVec_,&array );
          for( std::size_t i = 0; i < size(); i++ )
          {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%D %G\n",i,PetscRealPart(array[i]));
          }
          VecRestoreArray( ghostedVec_, &array );
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
          PetscSynchronizedFlush( PETSC_COMM_WORLD );
#else
          PetscSynchronizedFlush( PETSC_COMM_WORLD, PETSC_STDOUT );
#endif
      }

      // assign from other given PetscVector
      void assign( const ThisType& other )
      {
        // we want the 'other' to do all its communication right now before
        // we start copying values from it
        other.communicateIfNecessary();

        // copy vector entries
        ::Dune::Petsc::VecCopy( other.vec_, vec_ );

        // update ghost values
        updateGhostRegions();
      }

      // assign from other given SimpleBlockVector with same block size
      template <class Container>
      void assignVector( const SimpleBlockVector< Container, blockSize >& other )
      {
        Vec& vec = *getGhostedVector();

        // use mapping from the ghost mapper which removes duplicate dofs
        const auto& mapping = mappers_.ghostMapper().mapping();
        const size_t blocks = mapping.size();
        assert( blocks == other.size() );
        for( size_t b=0, bs = 0; b<blocks; ++b, bs += blockSize)
        {
          PetscInt block = mapping[ b ];
          const PetscScalar* values = static_cast< const PetscScalar* > (other.data()+bs);
          ::Dune::Petsc::VecSetValuesBlocked( vec, 1, &block, values, INSERT_VALUES );
        }
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );

        updateGhostRegions();
      }

      // assign from other given ISTLBlockVector with same block size
      template <class DofBlock>
      void assignVector( const ISTLBlockVector< DofBlock >& other )
      {
        assert( DofBlock :: dimension == blockSize );
        Vec& vec = *getGhostedVector();

        // use mapping from the ghost mapper which removes duplicate dofs
        const auto& mapping = mappers_.ghostMapper().mapping();
        const size_t blocks = mapping.size();
        assert( blocks == other.size() );
        for( size_t b=0; b<blocks; ++b )
        {
          PetscInt block = mapping[ b ];
          const PetscScalar* values = static_cast< const PetscScalar* > (&(other[ b ][ 0 ])) ;
          ::Dune::Petsc::VecSetValuesBlocked( vec, 1, &block, values, INSERT_VALUES );
        }
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );

        updateGhostRegions();
      }

      // assign from other given SimpleBlockVector with same block size
      template <class Container>
      void copyTo( SimpleBlockVector< Container, blockSize >& other ) const
      {
        typedef typename Container::FieldType FieldType;
        const PetscScalar *array = nullptr;
        VecGetArrayRead( ghostedVec_, &array );

        // use mapping from the ghost mapper which removes duplicate dofs
        const auto& mapping = mappers_.ghostMapper().mapping();
        const size_t blocks = mapping.size();
        assert( blocks == other.size() );
        for( size_t b=0; b<blocks; ++b )
        {
          const PetscScalar* petscBlock = array + (blockSize * mapping[ b ]);
          FieldType* otherBlock = other.data() + (b * blockSize);
          for( int i=0; i<blockSize; ++i )
          {
            otherBlock[ i ] = petscBlock[ i ];
          }
        }
        VecRestoreArrayRead( ghostedVec_, &array );
      }

      // assign from other given ISTLBlockVector with same block size
      template <class DofBlock>
      void copyTo ( ISTLBlockVector< DofBlock >& other ) const
      {
        assert( DofBlock :: dimension == blockSize );
        const PetscScalar *array = nullptr;
        VecGetArrayRead( ghostedVec_, &array );

        // use mapping from the ghost mapper which removes duplicate dofs
        const auto& mapping = mappers_.ghostMapper().mapping();
        const size_t blocks = mapping.size();
        assert( blocks == other.size() );
        for( size_t b=0; b<blocks; ++b )
        {
          const PetscScalar* petscBlock = array + (blockSize * mapping[ b ]);
          DofBlock& otherBlock = other[ b ];
          for( int i=0; i<blockSize; ++i )
          {
            otherBlock[ i ] = petscBlock[ i ];
          }
        }
        VecRestoreArrayRead( ghostedVec_, &array );
      }

      PetscVector& operator= ( const ThisType& other )
      {
        assign( other );
        return *this;
      }

      const MappersType &mappers() const { return mappers_; }

      std::size_t usedMemorySize() const
      {
        return size() * sizeof(value_type);
      }

      static void setMemoryFactor(const double memFactor)
      {
        // do nothing here
      }

      /** \brief move memory blocks backwards */
      void memMoveBackward( const int length, const int oldStartIdx, const int newStartIdx)
      {
        DUNE_THROW(NotImplemented,"memMoveBackward is to be implemented");
      }

      /** \brief move memory blocks forward */
      void memMoveForward(const int length, const int oldStartIdx, const int newStartIdx)
      {
        DUNE_THROW(NotImplemented,"memMoveForward is to be implemented");
      }

      void copyContent( const int newIndex, const int oldIndex )
      {
        DUNE_THROW(NotImplemented,"copyContent is to be implemented");
      }

    protected:
      // setup vector according to mapping sizes
      void init()
      {
        mappers_.update();

        const PetscInt localBlocks = mappers_.ghostMapper().interiorSize();
        const PetscInt numGhostBlocks = mappers_.ghostMapper().ghostSize();

        const PetscInt localSize = localBlocks * blockSize;
        const PetscInt globalSize = mappers_.parallelMapper().size() * blockSize;

        // finally, create the PETSc Vecs
        const PetscInt *ghostBlocks = mappers_.parallelMapper().mapping().data() + localBlocks;
        const auto& comm = mappers_.space().gridPart().comm();
        if( blockSize == 1 )
          ::Dune::Petsc::VecCreateGhost( comm, localSize, globalSize, numGhostBlocks, ghostBlocks, &vec_ );
        else
          ::Dune::Petsc::VecCreateGhostBlock( comm, blockSize, localSize, globalSize, numGhostBlocks, ghostBlocks, &vec_ );
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );
      }

      // delete vectors
      void removeObj()
      {
        ::Dune::Petsc::VecGhostRestoreLocalForm( vec_, &ghostedVec_ );
        ::Dune::Petsc::VecDestroy( &vec_ );
      }

      // communicate using the space default communication option
      void communicateIfNecessary () const
      {
        DefaultCommunicationOperationType op;
        communicateIfNecessary( op );
      }

      template <class Operation>
      void communicateIfNecessary (const Operation& op) const
      {
        // communicate this process' values
        if( communicateFlag_ && memorySequence_ < sequence_ )
        {
          if ( memorySequence_ < sequence_ )
          {
            if ( Operation::value == DFCommunicationOperation::add )
            {
              ::Dune::Petsc::VecGhostUpdateBegin( vec_, ADD_VALUES, SCATTER_REVERSE );
              ::Dune::Petsc::VecGhostUpdateEnd( vec_, ADD_VALUES, SCATTER_REVERSE );
            }
            ::Dune::Petsc::VecGhostUpdateBegin( vec_, INSERT_VALUES, SCATTER_FORWARD );
            ::Dune::Petsc::VecGhostUpdateEnd( vec_, INSERT_VALUES, SCATTER_FORWARD );

            memorySequence_ = sequence_;
          }

          communicateFlag_ = false;
        }
      }

      // Updates the ghost dofs, obtains them from the owning process
      void updateGhostRegions ()
      {
        ::Dune::Petsc::VecGhostUpdateBegin( vec_, INSERT_VALUES, SCATTER_FORWARD );
        ::Dune::Petsc::VecGhostUpdateEnd( vec_, INSERT_VALUES, SCATTER_FORWARD );
      }

      void vectorIsUpToDateNow () const
      {
        memorySequence_ = sequence_;
        communicateFlag_ = false;
      }

      /*
       * data fields
       */
      MappersType mappers_;
      Vec vec_;
      Vec ghostedVec_;

      mutable unsigned long memorySequence_ = 0;  // represents the state of the PETSc vec in memory
      mutable unsigned long sequence_ = 0;        // represents the modifications to the PETSc vec

      mutable bool communicateFlag_ = false;
      bool owner_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCVECTOR_HH
