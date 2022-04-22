#ifndef DUNE_FEM_DOFMANAGER_HH
#define DUNE_FEM_DOFMANAGER_HH

#include <cassert>
#include <string>
#include <list>

#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/version.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/interfaces.hh>
#endif

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/standardstreams.hh>
#include <dune/fem/misc/gridobjectstreams.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/datacollector.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/storage/dynamicarray.hh>
#include <dune/fem/storage/singletonlist.hh>

#include <dune/grid/common/datahandleif.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/ldbhandleif.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    /** @addtogroup DofManager

        @{
    **/

    // forward declaration
    template <class GridType> class DofManager;
    template <class DofManagerImp> class DofManagerFactory;


    ///////////////////////////////////////////////////////////////////////
    //
    //  ManagedIndexSetInterface
    //
    ///////////////////////////////////////////////////////////////////////
    /*! This class is the virtual interface for the index sets managed by
      the DofManager. The derived classes are of the type ManagedIndexSet<IndexSet>.
      This means we don't have to inherit every index set we want to use with
      this DofManager.
      */
    class ManagedIndexSetInterface
    {
      ManagedIndexSetInterface(const ManagedIndexSetInterface& org);
    protected:
      // use address of object as id
      typedef const void * IdentifierType;
      // pointer to compare index sets
      IdentifierType setPtr_;
      // reference counter
      size_t referenceCounter_;

      template< class IndexSet >
      explicit ManagedIndexSetInterface ( const IndexSet &iset )
      : setPtr_( getIdentifier( iset ) ),
        referenceCounter_( 1 )
      {
      }

    public:
      virtual ~ManagedIndexSetInterface () = default;

      //! resize of index set
      virtual void resize () = 0;
      //! compress of index set
      virtual bool compress () = 0;

      /** \copydoc Dune::PersistentObject :: backup */
      virtual void backup() const = 0;
      /** \copydoc Dune::PersistentObject :: restore */
      virtual void restore() = 0;

      //! new read/write methods using binary streams
      virtual void write( StandardOutStream& out ) const = 0;
      virtual void read( StandardInStream& out ) = 0;

      //! increase reference counter
      void addReference ( ) { ++referenceCounter_; }

      //! decrease reference counter and return true if zero reached
      bool removeReference ( )
      {
        return (--referenceCounter_ == 0);
      }

      template< class IndexSet >
      bool equals ( const IndexSet &iset ) const
      {
        return (getIdentifier( iset ) == setPtr_);
      }

    private:
      template< class IndexSet >
      IdentifierType getIdentifier ( const IndexSet &iset ) const
      {
        return static_cast< IdentifierType >( &iset );
      }
    };

    template <class IndexSetType, class EntityType> class RemoveIndicesFromSet;
    template <class IndexSetType, class EntityType> class InsertIndicesToSet;

    template <class IndexSetType, class EntityType>
    class ManagedIndexSet :
      public ManagedIndexSetInterface ,
      public LocalInlinePlus < ManagedIndexSet<IndexSetType,EntityType> , EntityType >
    {
      typedef LocalInterface<EntityType> LocalIndexSetObjectsType;

      static const bool isConsecutive =
        Capabilities::isConsecutiveIndexSet<IndexSetType>::v;

    protected:
      // the dof set stores number of dofs on entity for each codim
      IndexSetType & indexSet_;

      // insertion and removal of indices
      InsertIndicesToSet   <IndexSetType, EntityType> insertIdxObj_;
      RemoveIndicesFromSet <IndexSetType, EntityType> removeIdxObj_;

      LocalIndexSetObjectsType & insertList_;
      LocalIndexSetObjectsType & removeList_;

    public:
      //! type of base class
      typedef ManagedIndexSetInterface BaseType;

      //! Constructor of MemObject, only to call from DofManager
      ManagedIndexSet ( const IndexSetType & iset,
                        LocalIndexSetObjectsType & insertList,
                        LocalIndexSetObjectsType & removeList )
       : BaseType( iset )
       , indexSet_ (const_cast<IndexSetType &> (iset))
       , insertIdxObj_(indexSet_), removeIdxObj_(indexSet_)
       , insertList_(insertList)
       , removeList_(removeList)
      {
        this->setPtr_ = (void *) &indexSet_;

        if constexpr ( isConsecutive )
        {
          insertList_ += insertIdxObj_;
          removeList_ += removeIdxObj_;
        }
      }

      //! destructor
      ~ManagedIndexSet ()
      {
        if constexpr ( isConsecutive )
        {
          insertList_.remove( insertIdxObj_ );
          removeList_.remove( removeIdxObj_ );
        }
      }

      //! wrap resize of index set
      void resize ()
      {
        indexSet_.resize();
      }

      //! wrap compress of index set
      bool compress ()
      {
        return indexSet_.compress();
      }

      // forward backup call to indexSet
      virtual void backup () const
      {
        indexSet_.backup();
      }

      // forward restore call to indexSet
      virtual void restore ()
      {
        indexSet_.restore();
      }

      //! new write method
      virtual void read( StandardInStream& in ) { indexSet_.read( in ); }

      //! new write method
      virtual void write( StandardOutStream& out ) const { indexSet_.write( out ); }
    };

    /////////////////////////////////////////////////////////////
    //
    // DofStorageInterface
    //
    /////////////////////////////////////////////////////////////
    /** \brief Interface class for a dof storage object to be stored in
        discrete functions */
    class DofStorageInterface
    {
    protected:
      //! do not allow to create explicit instances
      DofStorageInterface() = default;

    public:
      typedef std::size_t SizeType;

      //! destructor
      virtual ~DofStorageInterface() = default;

      //! enable dof compression for dof storage (default is empty)
      virtual void enableDofCompression() { }

      //! size of space, i.e. mapper.size()
      virtual SizeType size () const = 0;
    };


    /** \brief Interface class for a dof storage object that can be managed
        (resized and compressed) by the DofManager
    */
    class ManagedDofStorageInterface : public DofStorageInterface
    {
    protected:
      //! do not allow to create explicit instances
      ManagedDofStorageInterface() = default;

    public:
      typedef typename DofStorageInterface::SizeType SizeType;
      // interface for MemObject lists
      typedef LocalInterface< SizeType > MemObjectCheckType;

      //! destructor
      virtual ~ManagedDofStorageInterface() = default;

      //! resize memory
      virtual void resize ( const bool enlargeOnly ) = 0;
      //! resize memory
      virtual void reserve (const SizeType newSize) = 0;
      //! compressed the underlying dof vector and clear the array if it's
      //! temporary and clearresizedarrays is enabled (default)
      virtual void dofCompress ( const bool clearResizedArrays ) = 0;
      //! return size of mem used by MemObject
      virtual size_t usedMemorySize() const = 0;
    };


    template <class MemObjectType> class ResizeMemoryObjects;
    template <class MemObjectType> class ReserveMemoryObjects;

    /*!
      A ManagedDofStorage holds the memory for one DiscreteFunction and the
      corresponding DofArrayMemory. ManagedDofStorageImplementation implements the
      basic features such as resize and dof compression. If a DiscreteFunction is signed in by a
      function space, then such a MemObject is created by the DofManager.
      The MemObject also knows the DofMapper from the function space which the
      discrete function belongs to. Here we dont know the exact type of the dof
      mapper therefore the methods newSize and calcInsertPoints of the mappers
      have to be virtual. This isnt a problem because this methods should only
      be called during memory reorganizing which is only once per timestep.
    */
    template <class GridImp, class MapperType , class DofArrayType>
    class ManagedDofStorageImplementation : public ManagedDofStorageInterface
    {
    public:
      typedef typename ManagedDofStorageInterface::SizeType SizeType;
    protected:
      typedef typename ManagedDofStorageInterface::MemObjectCheckType MemObjectCheckType;

      // type of this class
      typedef ManagedDofStorageImplementation <GridImp, MapperType , DofArrayType> ThisType;

      typedef DofManager<GridImp> DofManagerType;

      // reference to dof manager
      DofManagerType &dm_;

      // the dof set stores number of dofs on entity for each codim
      MapperType &mapper_;

      // Array which the dofs are stored in
      DofArrayType& array_;

      typedef ResizeMemoryObjects < ThisType > ResizeMemoryObjectType;
      typedef ReserveMemoryObjects  < ThisType > ReserveMemoryObjectType;
      ResizeMemoryObjectType  resizeMemObj_;
      ReserveMemoryObjectType reserveMemObj_;

      // true if data need to be compressed
      bool dataCompressionEnabled_;

    public:
      ManagedDofStorageImplementation(const ManagedDofStorageImplementation& ) = delete;

    protected:
      //! Constructor of ManagedDofStorageImplementation, only to call from derived classes
      ManagedDofStorageImplementation ( const GridImp& grid,
                                        const MapperType& mapper,
                                        DofArrayType& array )
        : dm_( DofManagerType :: instance( grid ) ),
          mapper_ ( const_cast<MapperType& >(mapper)),
          array_( array ),
          resizeMemObj_(*this),
          reserveMemObj_(*this),
          dataCompressionEnabled_(false)
      {
        // add to dof manager
        dm_.addDofStorage( *this );

        // set memory over estimate factor, only for DofArray
        array_.setMemoryFactor( dm_.memoryFactor() );
      }

      //! \brief destructor deleting MemObject from resize and reserve List
      ~ManagedDofStorageImplementation()
      {
        // remove from dof manager
        dm_.removeDofStorage( *this );
      }

    public:
      //! return object that calls resize of this memory object
      ResizeMemoryObjectType& resizeMemoryObject() { return resizeMemObj_; }

      //! return object that calls reserve of this memory object
      ReserveMemoryObjectType& reserveMemoryObject() { return reserveMemObj_; }

      //! return size of underlying array
      SizeType size () const override { return array_.size(); }

      //! resize the memory with the new size
      void resize ( const bool enlargeOnly ) override
      {
        resize( std::integral_constant< bool, Capabilities::isAdaptiveDofMapper< MapperType >::v >(), enlargeOnly );
      }

      //! reserve memory for what is coming
      inline void reserve ( const SizeType needed ) override
      {
        // if index set is compressible, then add requested size
        if( mapper().consecutive() )
        {
          const SizeType nSize = size() + (needed * SizeType(mapper().maxNumDofs()));
          array_.reserve( nSize );
        }
        else
        {
          // if compress is not needed just resize with given size
          // therefore use newSize to enlarge array
          assert( ! mapper().consecutive() );
          // resize array
          resize ( false );
        }
      }

      //! copy the dof from the rear section of the vector to the holes
      void dofCompress ( const bool clearResizedArrays ) override
      {
        // get current size
        const SizeType nSize = mapper().size();

        // if data is non-temporary do data compression
        if( dataCompressionEnabled_ )
        {
          // get consecutive information about mapper
          const bool consecutive = mapper().consecutive ();

          // get old size (which we still have in array)
          const SizeType oldSize = array_.size();

          // NOTE: new size can also be larger than old size
          // e.g. during loadBalancing when ghosts where
          // introduced before compressing the index set

          // begin with block zero since closing of holes
          // has to be done anyway if the mapper is consecutive
          const int numBlocks = mapper().numBlocks();
          for( int block = 0; block < numBlocks; ++block )
          {
            // move memory
            moveToFront( oldSize, block );

            // only close holes for consecutive mappers
            if( consecutive )
            {
              // run over all holes and copy array values to new place
              const SizeType holes = mapper().numberOfHoles( block );
              for( SizeType i = 0; i < holes; ++i )
              {
                const SizeType oldIndex = mapper().oldIndex( i, block );
                const SizeType newIndex = mapper().newIndex( i, block );

                assert( newIndex < nSize );
                // implements array_[ newIndex ] = array_[ oldIndex ] ;
                array_.copyContent( newIndex, oldIndex );
              }
            }
          }
        }

        // store new size, which should be smaller then actual size
        array_.resize( nSize );

        if( clearResizedArrays && ! dataCompressionEnabled_ )
        {
          // if enabled clear temporary data to avoid occasionally NaN values
          array_.clear();
        }
      }

      //! return used memory size
      size_t usedMemorySize() const override
      {
        return ((size_t) sizeof(ThisType) + array_.usedMemorySize());
      }

      //! enable dof compression for this MemObject
      void enableDofCompression() override
      {
        dataCompressionEnabled_ = true;
      }

      //! return reference to array for DiscreteFunction
      DofArrayType & getArray() { return array_; }

    protected:
      inline MapperType &mapper () const
      {
        return mapper_;
      }

      // resize for non-adaptive mappers
      void resize ( std::false_type, const bool enlargeOnly )
      {
        // note: The mapper might already have been updated, so do not use
        //       it to obtain old array size.
        mapper().update(); // make sure the mapper is up2date

        const SizeType newSize = mapper().size();
        const SizeType oldSize = array_.size();

        if( enlargeOnly && newSize < oldSize ) return ;

        if( newSize != oldSize )
          array_.resize( newSize );
      }

      // resize for adaptive mappers
      void resize ( std::true_type, const bool enlargeOnly )
      {
        // note: The mapper is adaptive and has been updated automatically, so
        //       do not use it to obtain old array size.
        const SizeType oldSize = array_.size();

        // get current size
        const SizeType nSize = mapper().size();

        // if enlarge only option is given only resize
        // if new size if larger than old size
        if( enlargeOnly && nSize <= oldSize ) return ;

        // if nothing changed do nothing
        if( nSize == oldSize ) return ;

        // resize memory to current value
        array_.resize( nSize );

        // if data is only temporary data, don't adjust memory
        if( ! dataCompressionEnabled_ || enlargeOnly ) return ;

        // now check all blocks beginning with the largest
        const int numBlocks = mapper().numBlocks();

        // initialize upperBound
        SizeType upperBound = oldSize ;

        // make sure offset of block 0 is zero
        assert( mapper().offSet( 0 ) == 0 );
        assert( mapper().oldOffSet( 0 ) == 0 );

        // skip block 0 (since offset == 0)
        for( int block = numBlocks-1; block >= 1; --block )
        {
          // get offsets
          const SizeType newOffSet = mapper().offSet( block );
          const SizeType oldOffSet = mapper().oldOffSet( block );

          // make sure new offset is larger
          assert( newOffSet >= oldOffSet );

          // if off set is not zero
          if( newOffSet > oldOffSet )
          {
            // calculate block size
            const SizeType blockSize = upperBound - oldOffSet;
            // move block backward
            array_.memMoveBackward( blockSize, oldOffSet, newOffSet );

            // update upper bound
            upperBound = oldOffSet;
          }
        }
      }

      // move array to rear insertion points
      void resizeAndMoveToRear ()
      {
      }

      //! move block to front again
      void moveToFront ( const SizeType oldSize, const int block )
      {
        // make sure offset of block 0 is zero
        assert( mapper().offSet( 0 ) == 0 );
        assert( mapper().oldOffSet( 0 ) == 0 );

        // block 0 has always offset 0
        if( block == 0 ) return;

        // get insertion point from block
        const SizeType oldOffSet = mapper().oldOffSet( block );

        // get new off set
        const SizeType newOffSet = mapper().offSet( block );

        // here we should have at least the same offsets
        assert( newOffSet <= oldOffSet );

        // only if block is not starting from zero
        if( newOffSet < oldOffSet )
        {
          // get number of blocks
          const int numBlocks = mapper().numBlocks();

          // for last section upperBound is size
          const SizeType upperBound
            = (block == numBlocks - 1) ? oldSize : mapper().oldOffSet( block + 1 );
          const SizeType blockSize = upperBound - oldOffSet;

          // move block forward
          array_.memMoveForward( blockSize, oldOffSet, newOffSet );
        }
      }
    };

    /*! A ManagedDofStorage holds the memory for one DiscreteFunction. */
    template <class GridImp, class MapperType , class DofArrayType>
    class ManagedDofStorage : public ManagedDofStorageImplementation< GridImp, MapperType, DofArrayType >
    {
      typedef ManagedDofStorageImplementation< GridImp, MapperType, DofArrayType > BaseType;
    protected:
      DofArrayType myArray_;
    public:
      //! Constructor of ManagedDofStorage
      ManagedDofStorage( const GridImp& grid,
                         const MapperType& mapper )
        : BaseType( grid, mapper, myArray_ ),
          myArray_( mapper.size() )
      {
      }
    };

    //! default implementation for creating a managed dof storage
    template< class DofStorageType, class GridType, class MapperType >
    static inline std::pair< DofStorageInterface* , DofStorageType* >
      allocateManagedDofStorage( const GridType& grid,
                                 const MapperType& mapper,
                                 const DofStorageType * = 0 )
    {
      // create managed dof storage
      typedef ManagedDofStorage< GridType, MapperType,
                                 DofStorageType > ManagedDofStorageType;

      ManagedDofStorageType* mds = new ManagedDofStorageType( grid, mapper );
      assert( mds );

      // return pair with dof storage pointer and array pointer
      return std::pair< DofStorageInterface* , DofStorageType* >
              ( mds , & mds->getArray () );
    }



    ///////////////////////////////////////////////////////////////
    //
    //  RestrictPorlong for Index Sets
    //
    ///////////////////////////////////////////////////////////////

    template <class IndexSetType, class EntityType>
    class RemoveIndicesFromSet
    : public LocalInlinePlus < RemoveIndicesFromSet<IndexSetType,EntityType> , EntityType >
    {
    private:
      // the dof set stores number of dofs on entity for each codim
      IndexSetType & indexSet_;

    public:
      // Constructor of MemObject, only to call from DofManager
      explicit RemoveIndicesFromSet ( IndexSetType & iset ) : indexSet_ (iset) {}

      //! apply wraps the removeEntity Method of the index set
      inline void apply ( EntityType & entity )
      {
        indexSet_.removeEntity( entity );
      }
    };

    template <class IndexSetType, class EntityType>
    class InsertIndicesToSet
    : public LocalInlinePlus < InsertIndicesToSet< IndexSetType, EntityType > , EntityType >
    {
    private:
      // the dof set stores number of dofs on entity for each codim
      IndexSetType & indexSet_;

    public:
      // Constructor of MemObject, only to call from DofManager
      explicit InsertIndicesToSet ( IndexSetType & iset ) : indexSet_ (iset) {}

      //! apply wraps the insertEntity method of the index set
      inline void apply ( EntityType & entity )
      {
        indexSet_.insertEntity( entity );
      }
    };

    template <class MemObjectType>
    class ResizeMemoryObjects
    : public LocalInlinePlus < ResizeMemoryObjects < MemObjectType > , typename MemObjectType::SizeType >
    {
    private:
      // the dof set stores number of dofs on entity for each codim
      MemObjectType& memobj_;

    public:
      typedef typename MemObjectType::SizeType SizeType;

      // Constructor of MemObject, only to call from DofManager
      ResizeMemoryObjects ( MemObjectType & mo ) : memobj_ (mo) {}
      ResizeMemoryObjects ( const ResizeMemoryObjects& org )
        : memobj_(org.memobj_)
      {}

      // resize mem object, parameter not needed
      inline void apply ( SizeType& enlargeOnly )
      {
        memobj_.resize( bool(enlargeOnly) );
      }
    };

    // this class is the object for a single MemObject to
    template <class MemObjectType>
    class ReserveMemoryObjects
    : public LocalInlinePlus < ReserveMemoryObjects < MemObjectType > , typename MemObjectType::SizeType >
    {
    private:
      // the dof set stores number of dofs on entity for each codim
      MemObjectType& memobj_;

    public:
      typedef typename MemObjectType::SizeType SizeType;

      // Constructor of MemObject, only to call from DofManager
      ReserveMemoryObjects ( MemObjectType & mo ) : memobj_ (mo) {}

      // reserve for at least chunkSize new values
      inline void apply ( SizeType& chunkSize )
      {
        memobj_.reserve( chunkSize );
      }
    };


    // this is the dofmanagers object which is being used during restriction
    // and prolongation process for adding and removing indices to and from
    // index sets which belong to functions that belong to that dofmanager
    template <class DofManagerType , class RestrictProlongIndexSetType, bool doResize >
    class IndexSetRestrictProlong  :
      public RestrictProlongInterface<
              RestrictProlongTraits<
                IndexSetRestrictProlong<DofManagerType,RestrictProlongIndexSetType,doResize>, double > >
    {
      DofManagerType & dm_;

      RestrictProlongIndexSetType & insert_;
      RestrictProlongIndexSetType & remove_;
    public:

      IndexSetRestrictProlong ( DofManagerType & dm , RestrictProlongIndexSetType & is, RestrictProlongIndexSetType & rm )
        : dm_(dm) , insert_( is ), remove_( rm ) {}

      //! restrict data to father and resize memory if doResize is true
      template <class EntityType>
      inline void restrictLocal ( const EntityType & father, const EntityType & son , bool initialize ) const
      {
        // insert index of father
        insert_.apply( father );
        // mark index of son for removal
        remove_.apply( son );

        // resize memory if doResize is true
        if ( doResize )
        {
          dm_.resizeMemory();
        }
      }

      template <class EntityType>
      inline void restrictFinalize( const EntityType &father ) const
      {}

      //! prolong data to children and resize memory if doResize is true
      template <class EntityType>
      inline void prolongLocal ( const EntityType & father, const EntityType & son , bool initialize ) const
      {
        // mark index of father for removal
        remove_.apply( father );
        // insert index of son
        insert_.apply( son );

        // resize memory if doResize is true
        if ( doResize )
        {
          dm_.resizeMemory();
        }
      }
    };

    // empty restrict prolong operator
    class EmptyIndexSetRestrictProlong  :
      public RestrictProlongInterface< RestrictProlongTraits< EmptyIndexSetRestrictProlong, double > >
    {
    public:
      EmptyIndexSetRestrictProlong() {}
      //! restrict data to father and resize memory if doResize is true
      template <class EntityType>
      inline void restrictLocal ( EntityType & father, EntityType & son , bool initialize ) const {}
      //! prolong data to children and resize memory if doResize is true
      template <class EntityType>
      inline void prolongLocal ( EntityType & father, EntityType & son , bool initialize ) const {}
    };


    class DofManError : public Exception {};

    /*!
     The DofManager is responsible for managing memory allocation and freeing
     for all discrete functions living on the grid the manager belongs to.
     There is only one DofManager per grid.
     Each discrete function knows its dofmanager and can sign in.
     If the grid is adapted, then the
     dofmanager reorganizes the memory if necessary. The DofManager holds a
     list of MemObjects which manage the memory and the corresponding
     mapper so they can determine the size of new memory.
     Furthermore the DofManager holds an IndexSet which the DofMapper needs for
     calculating the indices in the dof vector for a given entity and local dof
     number. This IndexSet is delivered to the mapper when a function space is
     created. The default value for the IndexSet is the DefaultIndexSet class
     which is mostly a wrapper for the grid indices.
    */
    // --DofManager
    template< class Grid >
    class DofManager :
#if HAVE_DUNE_ALUGRID
      public IsDofManager,
      public LoadBalanceHandleWithReserveAndCompress,
#endif
      // DofManager behaves like a communication data handle for load balancing
      public CommDataHandleIF< DofManager< Grid >, char >
    {
      typedef DofManager< Grid > ThisType;

      friend struct DefaultSingletonFactory< const Grid*, ThisType >;
      friend class DofManagerFactory< ThisType >;

    public:
      //! type of Grid this DofManager belongs to
      typedef Grid GridType;

      typedef typename GridObjectStreamTraits< GridType >::InStreamType  XtractStreamImplType;
      typedef typename GridObjectStreamTraits< GridType >::OutStreamType InlineStreamImplType;
    public:
      // types of inlining and xtraction stream types
      typedef MessageBufferIF< XtractStreamImplType > XtractStreamType;
      typedef MessageBufferIF< InlineStreamImplType > InlineStreamType;

      // types of data collectors
      typedef DataCollectorInterface<GridType, XtractStreamType >   DataXtractorType;
      typedef DataCollectorInterface<GridType, InlineStreamType  >  DataInlinerType;

      typedef typename GridType :: template Codim< 0 > :: Entity  ElementType ;

    private:
      typedef std::list< ManagedDofStorageInterface* > ListType;
      typedef typename ManagedDofStorageInterface::MemObjectCheckType MemObjectCheckType;
      typedef typename MemObjectCheckType::Traits::ParamType MemObjSizeType;
      typedef std::list< ManagedIndexSetInterface * > IndexListType;

      // list with MemObjects, for each DiscreteFunction we have one MemObject
      ListType memList_;

      // list of all different indexsets
      IndexListType indexList_;

      // the dofmanager belong to one grid only
      const GridType &grid_;

      // index set for mapping
      mutable DataInlinerType  dataInliner_;
      mutable DataXtractorType dataXtractor_;

      //! type of IndexSet change interfaces
      //// use const Entities as parameters (typedef here to avoid confusion)
      typedef const ElementType  ConstElementType;
      typedef LocalInterface< ConstElementType > LocalIndexSetObjectsType;

      mutable LocalIndexSetObjectsType insertIndices_;
      mutable LocalIndexSetObjectsType removeIndices_;

      // lists containing all MemObjects
      // to have fast access during resize and reserve
      mutable MemObjectCheckType resizeMemObjs_;
      mutable MemObjectCheckType reserveMemObjs_;

      //! if chunk size if small then defaultChunkSize is used
      const MemObjSizeType defaultChunkSize_;

      //! number of sequence, incremented every resize is called
      int sequence_;

    public:
      typedef IndexSetRestrictProlong< ThisType, LocalIndexSetObjectsType , true >
        NewIndexSetRestrictProlongType;
      typedef IndexSetRestrictProlong< ThisType, LocalIndexSetObjectsType , false >
        IndexSetRestrictProlongNoResizeType;

      // old type
      typedef EmptyIndexSetRestrictProlong IndexSetRestrictProlongType;

      // this class needs to call resizeMemory
      friend class IndexSetRestrictProlong< ThisType , LocalIndexSetObjectsType , true  > ;
      friend class IndexSetRestrictProlong< ThisType , LocalIndexSetObjectsType , false > ;

    private:
      // combine object holding all index set for restrict and prolong
      NewIndexSetRestrictProlongType indexSetRestrictProlong_;
      IndexSetRestrictProlongNoResizeType indexSetRestrictProlongNoResize_;

      // old type
      IndexSetRestrictProlongType indexRPop_;

      //! memory over estimation factor for re-allocation
      double memoryFactor_;

      //! true if temporary arrays should be reset to zero after resize
      const bool clearResizedArrays_;

      //**********************************************************
      //**********************************************************
      //! Constructor
      inline explicit DofManager ( const GridType *grid )
      : grid_( *grid ),
        defaultChunkSize_( 128 ),
        sequence_( 0 ),
        indexSetRestrictProlong_( *this, insertIndices_ , removeIndices_ ),
        indexSetRestrictProlongNoResize_( *this, insertIndices_ , removeIndices_ ),
        indexRPop_(),
        memoryFactor_( Parameter :: getValidValue( "fem.dofmanager.memoryfactor",  double( 1.1 ),
            [] ( double value ) { return value >= 1.0; } ) ),
        clearResizedArrays_( Parameter :: getValue("fem.dofmanager.clearresizedarrays", bool( true ) ) )
      {
        // only print memory factor if it deviates from the default value
        if( std::abs( memoryFactor_ - 1.1 ) > 1e-12 )
        {
          if( Parameter::verbose( Parameter::parameterOutput ) && (grid_.comm().rank() == 0) )
            std::cout << "Created DofManager with memory factor " << memoryFactor_ << "." << std::endl;
        }
      }

      //! Desctructor, removes all MemObjects and IndexSetObjects
      ~DofManager ();

    public:
      DofManager( const ThisType& ) = delete;

      //! return factor to over estimate new memory allocation
      double memoryFactor() const { return memoryFactor_; }

      /** \brief add index set to dof manager's list of index sets
       *
       *  During adaptation, all index sets known to the dof manager are notified
       *  of the changes.
       *
       *  To register an index set with the dof manager, it has to satisfy the
       *  following interface:
       *  \code
       *  void insertEntity ( const Element & );
       *  void removeEntity ( const Element & );
       *  void resize();
       *  bool compress();
       *  void write( OutStreamInterface<Traits>& );
       *  void read( InStreamInterface<Traits>& )
       *  \endcode
       *
       *  \param[in]  iset  index set to add to list
       */
      template <class IndexSetType>
      inline void addIndexSet (const IndexSetType &iset );

      /** \brief removed index set from dof manager's list of index sets
       *
       *  During adaptation, all index sets known to the dof manager are notified
       *  of the changes.
       *
       *  \param[in]  iset  index set to add to list
       */
      template <class IndexSetType>
      inline void removeIndexSet (const IndexSetType &iset );

      /** \brief add a managed dof storage to the dof manager.
          \param dofStorage  dof storage to add which must fulfill the
                 ManagedDofStorageInterface
      */
      template <class ManagedDofStorageImp>
      void addDofStorage(ManagedDofStorageImp& dofStorage);

      /** \brief remove a managed dof storage from the dof manager.
          \param dofStorage  dof storage to remove which must fulfill the
                 ManagedDofStorageInterface
      */
      template <class ManagedDofStorageImp>
      void removeDofStorage(ManagedDofStorageImp& dofStorage);

      //! returns the index set restriction and prolongation operator
      NewIndexSetRestrictProlongType & indexSetRestrictProlong ()
      {
        // hier muss statt dessen ein Combiniertes Object erzeugt werden.
        // dafuer sollte bei einhaengen der IndexSets ein Methoden Pointer
        // erzeugt werden, welcher die den IndexSet mit einem anderen Object
        // kombiniert
        return indexSetRestrictProlong_;
      }

      //! returns the index set restriction and prolongation operator
      IndexSetRestrictProlongNoResizeType& indexSetRestrictProlongNoResize()
      {
        // return index set restrict/prolong operator that is only inserting
        // and mark for removal indices but not doing resize
        return indexSetRestrictProlongNoResize_;
      }

      //! if dofmanagers list is not empty return true
      bool hasIndexSets() const
      {
        return ! insertIndices_.empty();
      }

      /** \brief return used memory size of all MemObjects in bytes. */
      size_t usedMemorySize () const
      {
        size_t used = 0;
        for(auto memObjectPtr : memList_)
          used += memObjectPtr->usedMemorySize();
        return used;
      }

      /** \brief resize memory before data restriction
          during grid adaptation is done.
      */
      void resizeForRestrict ()
      {
        resizeMemory();
      }

      /** \brief reserve memory for at least nsize elements,
       *         dummy is needed for dune-grid ALUGrid version */
      void reserveMemory ( std::size_t nsize, bool dummy = false )
      {
        MemObjSizeType localChunkSize =
          std::max( MemObjSizeType(nsize), defaultChunkSize_ );
        assert( localChunkSize > 0 );

        // reserves (size + chunkSize * elementMemory), see above
        reserveMemObjs_.apply ( localChunkSize );
      }

      /** \brief return number of sequence, if dofmanagers memory was changed by
          calling some method like resize, then also this number will increase
         \note The increase of this number could be larger than 1
      */
      int sequence () const { return sequence_; }

      /** \brief Resize index sets and memory due to what the mapper has as new size.
          \note This will increase the sequence counter by 1.
      */
      void resize()
      {
        for(auto indexSetPtr : indexList_)
          indexSetPtr->resize();
        resizeMemory();
      }

      /** \brief Inserts entity to all index sets added to dof manager. */
      inline void insertEntity( ConstElementType & element )
      {
        // insert new index
        insertIndices_.apply( element );

        // resize memory
        resizeMemory();
      }

      /** \brief Removes entity from all index sets added to dof manager. */
      inline void removeEntity( ConstElementType & element )
      {
        removeIndices_.apply( element );
      }

      //! resize the MemObject if necessary
      void resizeMemory()
      {
        MemObjSizeType enlargeOnly( 0 );
        // pass dummy parameter
        resizeMemObjs_.apply ( enlargeOnly );
      }

      //! resize the MemObject if necessary
      void enlargeMemory()
      {
        MemObjSizeType enlargeOnly( 1 );
        // pass dummy parameter
        resizeMemObjs_.apply ( enlargeOnly );
      }

      /** \brief increase the DofManagers internal sequence number
          \note  This will increase the sequence counter by 1.
      */
      void incrementSequenceNumber ()
      {
        // mark next sequence
        ++sequence_;

        // check that sequence number is the same for all processes
        assert( sequence_ == grid_.comm().max( sequence_ ) );
      }

      //- --compress
      /** \brief Compress all data that is hold by this dofmanager
          \note  This will increase the sequence counter by 1.
      */
      void compress()
      {
        // mark next sequence
        incrementSequenceNumber ();

        // compress indexsets first
        for(auto indexSetPtr : indexList_)
        {
          // reset compressed so the next time compress of index set is called
          indexSetPtr->compress();
        }

        // compress all data now
        for(auto memObjectPtr : memList_)
        {
          // if correponding index was not compressed yet, this is called in
          // the MemObject dofCompress, if index has not changes, nothing happens
          // if IndexSet actual needs  no compress, nothing happens to the
          // data either
          // also data is resized, which means the vector is getting shorter
          memObjectPtr->dofCompress ( clearResizedArrays_ );
        }
      }

      //! communicate new sequence number
      bool notifyGlobalChange( const bool wasChanged ) const
      {
        // make sure that wasChanged is the same on all cores
        int wasChangedCounter = int( wasChanged );
        return bool( grid_.comm().max( wasChangedCounter ) );
      }

      //! add data handler for data inlining to dof manager
      template <class DataCollType>
      void addDataInliner ( DataCollType & d)
      {
        dataInliner_ += d;
      }

      //! clear data inliner list
      void clearDataInliners ()
      {
        dataInliner_.clear();
      }

      //! add data handler for data xtracting to dof manager
      template <class DataCollType>
      void addDataXtractor ( DataCollType & d)
      {
        dataXtractor_ += d;
      }

      //! clear data xtractor list
      void clearDataXtractors ()
      {
        dataXtractor_.clear();
      }

      //////////////////////////////////////////////////////////
      //  CommDataHandleIF methods
      //////////////////////////////////////////////////////////

      //! the dof manager only transfers element data during load balancing
      bool contains( const int dim, const int codim ) const
      {
        return ( codim == 0 );
      }

      //! fixed size is false
      bool fixedSize( const int dim, const int codim ) const
      {
        return false;
      }

      //! for convenience
      template <class Entity>
      size_t size( const Entity& ) const
      {
        DUNE_THROW(NotImplemented,"DofManager::size should not be called!");
        return 0;
      }

      //! packs all data attached to this entity
      void gather( InlineStreamType& str, ConstElementType& element ) const
      {
        dataInliner_.apply(str, element);

        // remove entity from index sets
        const_cast< ThisType & >( *this ).removeEntity( element );
      }

      template <class MessageBuffer, class Entity>
      void gather( MessageBuffer& str, const Entity& entity ) const
      {
        DUNE_THROW(NotImplemented,"DofManager::gather( entity ) with codim > 0 not implemented");
      }

      //! unpacks all data attached of this entity from message buffer
      void scatter ( XtractStreamType& str, ConstElementType& element, size_t )
      {
        // insert entity into index sets
        insertEntity( element );

        // here the elements already have been created
        // that means we can xtract data
        dataXtractor_.apply(str, element);
      }

      //! unpacks all data of this entity from message buffer
      template <class MessageBuffer, class Entity>
      void scatter ( MessageBuffer & str, const Entity& entity, size_t )
      {
        DUNE_THROW(NotImplemented,"DofManager::scatter( entity ) with codim > 0 not implemented");
      }

      //********************************************************
      // Interface for PersistentObjects
      //********************************************************

      /** \copydoc Dune::PersistentObject :: backup */
      void backup () const
      {
        for(auto indexSetPtr : indexList_)
          indexSetPtr->backup();
      }

      /** \copydoc Dune::PersistentObject :: restore */
      void restore ()
      {
        for(auto indexSetPtr : indexList_)
          indexSetPtr->restore();

        // make all index sets consistent
        // before any data is read this can be
        // assured since DofManager is an
        // AutoPersistentObject and thus in the
        // beginning of the list, fater the grid of course
        resize();
      }

      //********************************************************
      // read/write using fem streams
      //********************************************************
      /** \brief write all index sets to a given stream
       *
       *  \param[out]  out  stream to write to
       */
      template < class OutStream >
      void write( OutStream& out ) const
      {
        for(auto indexSetPtr : indexList_)
          indexSetPtr->write( out );
      }

      /** \brief read all index sets from a given stream
       *
       *  \param[in] in  stream to read from
       */
      template < class InStream >
      void read( InStream& in )
      {
        for(auto indexSetPtr : indexList_)
          indexSetPtr->read( in );
      }

      //********************************************************
      // Interface for DofManager access
      //********************************************************

      /** \brief obtain a reference to the DofManager for a given grid
       *
       *  \param[in]  grid  grid for which the DofManager is desired
       *
       *  \returns a reference to the singleton instance of the DofManager
       */
      static inline ThisType& instance( const GridType& grid )
      {
        typedef DofManagerFactory< ThisType > DofManagerFactoryType;
        return DofManagerFactoryType :: instance( grid );
      }
    };

    //***************************************************************************
    //
    //  inline implemenations
    //
    //***************************************************************************

    template <class GridType>
    inline DofManager<GridType>::~DofManager ()
    {
      // enable output if verbosity level is debugOutput
      const bool verbose = Parameter::verbose( Parameter::debugOutput );
      if(memList_.size() > 0)
      {
        while( memList_.rbegin() != memList_.rend())
        {
          DofStorageInterface * mobj = (* memList_.rbegin() );
          if( verbose )
            std::cout << "Removing '" << mobj << "' from DofManager!\n";
          memList_.pop_back();
        }
      }

      if(indexList_.size() > 0)
      {
        if( verbose )
          std::cerr << "ERROR: Not all index sets have been removed from DofManager yet!" << std::endl;
        while ( indexList_.rbegin() != indexList_.rend())
        {
          ManagedIndexSetInterface* iobj = (* indexList_.rbegin() );
          indexList_.pop_back();
          if(iobj) delete iobj;
        }
      }
    }

    template <class GridType>
    template <class IndexSetType>
    inline void DofManager<GridType>::
    addIndexSet (const IndexSetType &iset )
    {
      // only call in single thread mode
      if( ! Fem :: MPIManager :: singleThreadMode() )
      {
        assert( Fem :: MPIManager :: singleThreadMode() );
        DUNE_THROW(InvalidStateException,"DofManager::addIndexSet: only call in single thread mode!");
      }

      typedef ManagedIndexSet< IndexSetType, ConstElementType > ManagedIndexSetType;
      ManagedIndexSetType * indexSet = 0;

      // search index set list in reverse order to find latest index sets faster
      auto endit = indexList_.rend();
      for(auto it = indexList_.rbegin(); it != endit; ++it )
      {
        ManagedIndexSetInterface *set = *it;
        if( set->equals( iset ) )
        {
          set->addReference();

          indexSet = static_cast< ManagedIndexSetType * >( set );
          break;
        }
      }

      if( !indexSet )
      {
        indexSet = new ManagedIndexSetType ( iset, insertIndices_ , removeIndices_  );
        indexList_.push_back( static_cast< ManagedIndexSetInterface * >( indexSet ) );
      }
    }

    template <class GridType>
    template <class IndexSetType>
    inline void DofManager<GridType>::removeIndexSet ( const IndexSetType &iset )
    {
      // only call in single thread mode
      if( ! Fem :: MPIManager :: singleThreadMode() )
      {
        assert( Fem :: MPIManager :: singleThreadMode() );
        DUNE_THROW(InvalidStateException,"DofManager::removeIndexSet: only call in single thread mode!");
      }

      // search index set list in reverse order to find latest index sets faster
      auto endit = indexList_.rend();
      for( auto it = indexList_.rbegin(); it != endit; ++it )
      {
        ManagedIndexSetInterface *set = *it;
        if( set->equals( iset ) )
        {
          if( set->removeReference() )
          {
            // reverse iterators cannot be erased directly, so erase the base
            // (forward) iterator
            // Note: see, e.g., Stroustrup, section 16.3.2 about the decrement
            auto fit = it.base();
            indexList_.erase( --fit );
            // delete proxy
            delete set;
          }
          return;
        }
      }

      // we should never get here
      DUNE_THROW(InvalidStateException,"Could not remove index from DofManager set!");
    }

    template <class GridType>
    template <class ManagedDofStorageImp>
    void DofManager<GridType>::addDofStorage(ManagedDofStorageImp& dofStorage)
    {
      // make sure we got an ManagedDofStorage
      ManagedDofStorageInterface* obj = &dofStorage;

      // push_front, makes search faster
      memList_.push_front( obj );

      // add the special object to the memResize list object
      resizeMemObjs_  += dofStorage.resizeMemoryObject();

      // the same for the reserve call
      reserveMemObjs_ += dofStorage.reserveMemoryObject();
    }


    template <class GridType>
    template <class ManagedDofStorageImp>
    void DofManager<GridType>::removeDofStorage(ManagedDofStorageImp& dofStorage)
    {
      // make sure we got an ManagedDofStorage
      auto obj = &dofStorage;

      // search list starting from tail
      auto endit = memList_.end();
      for( auto it = memList_.begin();it != endit ; ++it)
      {
        if(*it == obj)
        {
          // alloc new mem and copy old mem
          memList_.erase( it );

          // remove from list
          resizeMemObjs_.remove( dofStorage.resizeMemoryObject() );
          reserveMemObjs_.remove( dofStorage.reserveMemoryObject() );

          return ;
        }
      }
    }

    //@}


    /** \class DofManagerFactory
     *  \ingroup DofManager
     *  \brief Singleton provider for the DofManager
     *
     *  DofManagerFactory guarantees that at most one instance of DofManager
     *  is generated for each grid.
     */
    template< class DofManagerImp >
    class DofManagerFactory
    {
      typedef DofManagerFactory< DofManagerImp > ThisType;

    public:
      typedef DofManagerImp DofManagerType;
      typedef typename DofManagerType :: GridType GridType;

    private:
      typedef const GridType *KeyType;

      typedef SingletonList< KeyType, DofManagerType > DMProviderType;

      // declare friendship becase of methods instance
      friend class DofManager< GridType >;

    protected:
      /** \brief obtain a reference to the DofManager for a given grid
       *
       *  \param[in]  grid  grid for which the DofManager is desired
       *
       *  \returns a reference to the singleton instance of the DofManager
       */
      inline static DofManagerType &instance ( const GridType &grid )
      {
        DofManagerType *dm = getDmFromList( grid );
        if( !dm )
          return DMProviderType :: getObject( &grid );
        return *dm;
      }

      //! writes DofManager of corresponding grid, when DofManager exists
      inline static bool
      writeDofManagerNew ( const GridType &grid,
                           const std :: string &filename,
                           int timestep )
      {
        DofManagerType *dm = getDmFromList( grid );
        /*
        if( dm )
          return dm->writeIndexSets( filename, timestep );
          */
        return false;
      }

      //! reads DofManager of corresponding grid, when DofManager exists
      inline static bool
      readDofManagerNew ( const GridType &grid,
                          const std :: string &filename,
                          int timestep )
      {
        DofManagerType *dm = getDmFromList( grid );
        /*
        if( dm )
          return dm->readIndexSets( filename, timestep );
          */
        return false;
      }

    public:
      //! delete the dof manager that belong to the given grid
      inline static void deleteDofManager ( DofManagerType &dm )
      {
        DMProviderType :: removeObject( &dm );
      }

    private:
      // return pointer to dof manager for given grid
      inline static DofManagerType *getDmFromList( const GridType &grid )
      {
        return (DMProviderType :: getObjFromList( &grid )).first;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMANAGER_HH
