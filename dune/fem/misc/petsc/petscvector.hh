// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCVECTOR_HH
#define DUNE_FEM_PETSCVECTOR_HH

#include <dune/fem/storage/envelope.hh> 

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscslavedofprovider.hh>
#include <dune/fem/misc/petsc/petscghostarraybuilder.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscdofblock.hh>
#include <dune/common/static_assert.hh>


namespace Dune 
{

  namespace Fem 
  {

    // forward declarations
    template< typename >      class PetscDofBlock;
    template< class DFSpace > class PetscVector;

    /*! ManagedDofStorage for PetscDiscreteFunction using PetscVector */
    template < class DiscreteFunctionSpace, class Mapper >
    class PetscManagedDofStorage 
      : public ManagedDofStorageImplementation< typename DiscreteFunctionSpace :: GridType, 
                                                Mapper,
                                                PetscVector< DiscreteFunctionSpace > > 
    {
      typedef typename DiscreteFunctionSpace :: GridType GridType;
      typedef PetscSlaveDofProvider< DiscreteFunctionSpace > PetscSlaveDofProviderType;
      typedef Mapper MapperType ;
      typedef PetscVector< DiscreteFunctionSpace > DofArrayType ;
      typedef ManagedDofStorageImplementation< GridType, MapperType, DofArrayType >  BaseType;
    protected:
      DofArrayType myArray_;
    public:
      //! Constructor of ManagedDofStorageImpl, only to call from DofManager 
      PetscManagedDofStorage( const DiscreteFunctionSpace& space, 
                              const MapperType& mapper, 
                              const std::string& name )
        : BaseType( space.grid(), mapper, name, myArray_ ),
          myArray_( space )
      {
      }
    };


    /*! specialization of SpecialArrayFeatures for PetscVector 
     * dealing with the strange PetscVec */
    template< class DFS > 
    struct SpecialArrayFeatures< PetscVector< DFS > > 
    {
      typedef PetscVector< DFS >  ArrayType ;
      /** \brief value type of array, i.e. double */
      typedef typename ArrayType :: value_type ValueType;

      /** \brief return used memory size of Array */
      static size_t used(const ArrayType & array)
      {
        return array.size() * sizeof(ValueType);
      }

      /** \brief set memory overestimate factor, here does nothing */
      static void setMemoryFactor(ArrayType & array, const double memFactor)
      {
        // do nothing here
      }

      /** \brief move memory blocks backwards */
      static void memMoveBackward(ArrayType& array, const int length,
                                  const int oldStartIdx, const int newStartIdx)
      {
        DUNE_THROW(NotImplemented,"memMoveBackward is to be implemented");
      }

      /** \brief move memory blocks forward */
      static void memMoveForward(ArrayType& array, const int length,
                                 const int oldStartIdx, const int newStartIdx)
      {
        DUNE_THROW(NotImplemented,"memMoveForward is to be implemented");
      }

      static void assign( ArrayType& array, const int newIndex, const int oldIndex )
      {
        /*
        typedef typename ArrayType :: DofBlockPtrType DofBlockPtrType;
        DofBlockPtrType newBlock = array.block( newIndex );
        DofBlockPtrType oldBlock = array.block( oldIndex );

        const unsigned int blockSize = ArrayType :: blockSize;
        for( unsigned int i = 0; i < blockSize; ++i )
          (*newBlock)[ i ] = (*oldBlock)[ i ];
          */
      }
    };


    /* ========================================
     * class PetscVector
     * =======================================
     */
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

    public:
      typedef PetscSlaveDofProvider< DFSpace > PetscSlaveDofsType;
      typedef PetscScalar  value_type ;

      static const int blockSize = DFSpace :: localBlockSize;
      typedef typename PetscSlaveDofsType :: PetscDofMappingType  PetscDofMappingType;
      
      typedef PetscDofBlock< ThisType >                       DofBlockType;
      typedef PetscDofBlock< const ThisType >                 ConstDofBlockType;
      typedef typename DofBlockType::DofIterator              DofIteratorType;
      typedef typename ConstDofBlockType::DofIterator         ConstDofIteratorType;
      typedef Envelope< DofBlockType >                        DofBlockPtrType; 
      typedef Envelope< ConstDofBlockType >                   ConstDofBlockPtrType;
      typedef typename DofBlockType::IndexType                IndexType;

      typedef typename DFSpace::template CommDataHandle<void>::OperationType CommunicationOperationType;

      PetscVector ( const DFSpace& dfSpace )
      : petscSlaveDofs_( dfSpace ),
        memorySequence_( 0 ),
        sequence_( 0 ),
        communicateFlag_( false ),
        localSize_( 0 ),
        numGhosts_( 0 )
      {
        dune_static_assert( CommunicationOperationType::value == DFCommunicationOperation::copy ||
                            CommunicationOperationType::value == DFCommunicationOperation::add,
                            "only copy/add are available communication operations for petsc");
        // init vector 
        init();
      }

      // TODO: think about sequence overflows...
      PetscVector ( const ThisType &other )
      : petscSlaveDofs_( other.petscSlaveDofs_.space() ),
        memorySequence_( 0 ),
        sequence_( 0 ),
        communicateFlag_( false ),
        localSize_( other.localSize_ ),
        numGhosts_( other.numGhosts_ )
      {
        // assign vectors 
        assign( other );
      }
      
      ~PetscVector ()
      {
        // destroy vectors 
        removeObj();
      }

      size_t size () const { return static_cast< size_t >( localSize_ + numGhosts_ ); }

      void resize( const size_t newsize ) 
      {
        /*
        std::vector< double > values( dofMapping().size() );
        typedef typename std::vector< double > :: iterator  iterator ;

        const DofIteratorType end = dend ();
        iterator value = values.begin();
        for( DofIteratorType it = dbegin(); it != end ; ++ it, ++value ) 
        {
          if( value == values.end() ) break ;
          assert( value != values.end() );
          *value = *it ;
        }
        */

        // TODO: keep old data stored in current vector
        // remove old vectors 
        removeObj();

        // initialize new 
        init();

        /*
        const size_t vsize = std::min( values.size(), size() );
        DofIteratorType it = dbegin();
        for( size_t i=0; i<vsize; ++ i, ++ it ) 
          *it = values[ i ];      
        */

        hasBeenModified ();
      }

      void reserve( const size_t capacity ) 
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

      // force communication _now_
      void communicateNow () const
      {
        communicateFlag_ = true;
        ++sequence_;
        communicateIfNecessary();
      }

      //DofBlockPtrType operator [] ( const IndexType index ) { return block( index ); }
      //ConstDofBlockPtrType operator [] const ( const IndexType index ) { return block( index ); }

      DofBlockPtrType block ( IndexType index ) 
      {
        assert( index < dofMapping().size() );
        return DofBlockPtrType( typename DofBlockType::UnaryConstructorParamType( *this, index ) );
      }

      ConstDofBlockPtrType block ( IndexType index ) const
      {
        assert( index < dofMapping().size() );
        return ConstDofBlockPtrType( typename ConstDofBlockType::UnaryConstructorParamType( *this, index ) );
      }

      DofIteratorType dbegin () { return DofIteratorType( *this, 0, 0 ); }
      ConstDofIteratorType dbegin () const { return ConstDofIteratorType( *this, 0, 0 ); }
      DofIteratorType dend() { return DofIteratorType( *this, dofMapping().size() ); }
      ConstDofIteratorType dend() const { return ConstDofIteratorType( *this, dofMapping().size() ); }

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
          for( int i=0; i < localSize_ + numGhosts_; i++ ) 
          {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%D %G\n",i,PetscRealPart(array[i]));
          }
          VecRestoreArray( ghostedVec_, &array );
          PetscSynchronizedFlush( PETSC_COMM_WORLD );
      }

    protected:
      // setup vector according to mapping sizes 
      void init() 
      {
        // set up the DofMapping instance and all variables depending on it
        localSize_ = dofMapping().numOwnedDofBlocks() * blockSize;
        numGhosts_ = dofMapping().numSlaveBlocks()    * blockSize;

        //std::cout << "PetscVector::init: "<< localSize_ << "  " << numGhosts_ << std::endl;
        assert( static_cast< size_t >( localSize_ + numGhosts_ ) == dofMapping().size() * blockSize );

        // set up the ghost array builder
        typedef PetscGhostArrayBuilder< PetscSlaveDofsType, PetscDofMappingType > PetscGhostArrayBuilderType;
        PetscGhostArrayBuilderType ghostArrayBuilder( petscSlaveDofs_, dofMapping() );
        assert( int( ghostArrayBuilder.size() ) == dofMapping().numSlaveBlocks() );
        
        // finally, create the PETSc Vecs
        ::Dune::Petsc::VecCreateGhost( 
          static_cast< PetscInt >( localSize_ ), 
          PETSC_DECIDE, 
          static_cast< PetscInt >( numGhosts_ ), 
          ghostArrayBuilder.array(),
          &vec_ 
        );
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );
      }

      // delete vectors 
      void removeObj() 
      {
        ::Dune::Petsc::VecGhostRestoreLocalForm( vec_, &ghostedVec_ );
        ::Dune::Petsc::VecDestroy( &vec_ );
      }

      // assign from other given PetscVector
      void assign( const ThisType& other ) 
      {
        // we want the 'other' to do all its communication right now before
        // we start copying values from it
        other.communicateIfNecessary();

        // Do the copying on the PETSc level
        ::Dune::Petsc::VecDuplicate( other.vec_, &vec_ );
        ::Dune::Petsc::VecCopy( other.vec_, vec_ );
        ::Dune::Petsc::VecGhostGetLocalForm( vec_, &ghostedVec_ );

        updateGhostRegions();
      }

      PetscVector ();
      PetscVector& operator= ( const ThisType& );

      PetscDofMappingType& dofMapping () { return petscSlaveDofs_.dofMapping(); }
      const PetscDofMappingType& dofMapping () const { return petscSlaveDofs_.dofMapping(); }
      
      void communicateIfNecessary () const
      {
        // communicate this process' values
        if( communicateFlag_ && memorySequence_ < sequence_ )
        {
          if ( memorySequence_ < sequence_ )
          {
            if ( CommunicationOperationType::value == DFCommunicationOperation::add )
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
      PetscSlaveDofsType petscSlaveDofs_;
      Vec vec_;
      Vec ghostedVec_;

      mutable unsigned long memorySequence_; // represents the state of the PETSc vec in memory
      mutable unsigned long sequence_; // represents the the modifications to the PETSc vec

      mutable bool communicateFlag_;
      PetscInt localSize_;
      PetscInt numGhosts_;

    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC 

#endif // DUNE_FEM_PETSCVECTOR_HH
