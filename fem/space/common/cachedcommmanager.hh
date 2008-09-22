#ifndef DUNE_CACHED_COMMUNICATION_MANAGER_HH
#define DUNE_CACHED_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <queue>
#include <vector>

//- Dune includes  
#include <dune/common/misc.hh>
#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_ALUGRID
// inlcude alugrid to have to communicator class from ALUGrid 
#include <dune/grid/alugrid.hh>
#endif

//- Dune-fem includes 
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/entitycommhelper.hh>
#include <dune/fem/space/common/commindexmap.hh>

namespace Dune
{
  
/** @addtogroup Communication Communication 
    @{
**/
  
  // only if ALUGrid found and was build for parallel runs 
#if HAVE_ALUGRID && ALU3DGRID_PARALLEL 
  //! class to build up a map of all dofs of entities to be exchanged
  //! during a communication procedure. This speeds up the communication
  //! procedure, because no grid traversal is necessary to exchange data.
  //! this class is singleton for different discrete function space,
  //! because the dof mapping is always the same.
  template< class Space >
  class DependencyCache 
  {
  public:
    //! type of discrete function space 
    typedef Space SpaceType;

     //! type of grid part 
    typedef typename SpaceType :: GridPartType GridPartType; 
   
    class UnsaveObjectStream;

  protected:
    // type of communication indices 
    typedef CommunicationIndexMap IndexMapType;

    // type of IndexMapVector
    typedef IndexMapType *IndexMapVectorType;

    // type of set of links 
    typedef std :: set< int >  LinkStorageType;

    // ALUGrid send/recv buffers
    typedef ALU3DSPACE ObjectStream ObjectStreamType;
    
    // type of communicator
    typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType;
    // type of communication implementation
    typedef ALU3DSPACE MpAccessMPI MPAccessImplType;
    
     //! type of communication buffer vector 
    typedef std :: vector< ObjectStreamType > ObjectStreamVectorType;
  
  protected:
    const SpaceType &space_; 
    const GridPartType &gridPart_;

    const InterfaceType interface_;
    const CommunicationDirection dir_; 
    
    const int myRank_;
    const int mySize_; 
    
    LinkStorageType linkStorage_; 

    IndexMapType *recvIndexMap_;
    IndexMapType *sendIndexMap_;

    // vector containing the links of this process 
    std :: vector< int > linkRank_;

    // ALUGrid communicatior Class 
    MPAccessInterfaceType *mpAccess_;

    //! communication buffers 
    ObjectStreamVectorType buffer_;

    //! number of links 
    int nLinks_;

    //! know grid sequence number 
    int sequence_; 
    
    // exchange time 
    double exchangeTime_;
    // setup time 
    double buildTime_;

  protected:
    template< class LinkStorage, class IndexMapVector, InterfaceType CommInterface >
    class LinkBuilder;

  public:
    //! constructor taking space 
    DependencyCache ( const SpaceType &space,
                      const InterfaceType interface, 
                      const CommunicationDirection dir) 
    : space_( space ),
      gridPart_( space_.gridPart() ),
      interface_( interface ),
      dir_(dir),
      myRank_( gridPart_.grid().comm().rank() ),
      mySize_( gridPart_.grid().comm().size() ),
      linkStorage_(),
      recvIndexMap_( new IndexMapType[ mySize_ ] ),
      sendIndexMap_( new IndexMapType[ mySize_ ] ),
      linkRank_(),
      // create mpAccess with world communicator  
      // only when size > 1 
      // mpAccess_( (mySize_ > 1) ? 
      //    (new MPAccessImplType( MPIHelper::getCommunicator() )) : 0),
      mpAccess_( new MPAccessImplType( MPIHelper::getCommunicator() ) ),
      nLinks_( 0 ),
      sequence_( -1 ),
      exchangeTime_( 0.0 ),
      buildTime_( 0.0 )
    {
    }

  private:
    // prohibit copying
    DependencyCache ( const DependencyCache & );

  public:
    //! destrcutor removeing mpAccess 
    inline ~DependencyCache ()
    {
      delete mpAccess_;
      mpAccess_ = 0;

      delete [] sendIndexMap_;
      sendIndexMap_ = 0;
      
      delete [] recvIndexMap_;
      recvIndexMap_ = 0;
    }

  public:
    //! return communication interface 
    InterfaceType communicationInterface() const {
      return interface_;
    }

    //! return communcation direction 
    CommunicationDirection communicationDirection() const
    {
      return dir_;
    }

    //! return time needed for last build 
    double buildTime() const { return buildTime_; }

    //! return time needed for last exchange  
    double exchangeTime() const { return exchangeTime_; }

    // build linkage and index maps 
    inline void buildMaps ();

  protected:
    // check consistency of maps 
    inline void checkConsistency ();

    template< class LS, class IMV, InterfaceType CI >
    inline void buildMaps ( LinkBuilder< LS, IMV, CI > &handle );

  public:
    //! return number of links
    inline int nlinks () const
    {
      return nLinks_;
    }

    //! check if grid has changed and rebuild cache if necessary 
    inline void rebuild () 
    {
      // only in parallel we have to do something 
      if( mySize_ <= 1 ) return;

#ifndef NDEBUG
      // make sure buildMaps is called on every process 
      // otherwise the programs wait here until forever 
      int willRebuild = (sequence_ != space_.sequence()) ? 1 : 0;
      const int myRebuild = willRebuild;

      // send willRebuild from rank 0 to all 
      gridPart_.grid().comm().broadcast( &willRebuild, 1 , 0);

      assert( willRebuild == myRebuild );
#endif

      // check whether grid has changed.
      if( sequence_ != space_.sequence() )
      {
        // take timer needed for rebuild 
        Timer buildTime;

        buildMaps();
        sequence_ = space_.sequence();

        // store time needed 
        buildTime_ = buildTime.elapsed();
      }
    }
      
    //! exchange data of discrete function 
    template< class DiscreteFunction, class Operation >
    inline void exchange ( DiscreteFunction &discreteFunction,
                           const Operation *operation );

    //! write data of discrete function to buffer 
    template< class DiscreteFunction >
    inline void writeBuffer ( ObjectStreamVectorType &osv,
                              const DiscreteFunction &discreteFunction ) const;
    
    //! read data of discrete function from buffer  
    template< class DiscreteFunctionType, class Operation >
    inline void readBuffer ( ObjectStreamVectorType &osv,
                             DiscreteFunctionType &discreteFunction,
                             const Operation *operation ) const;
    
    //! return reference to mpAccess object
    inline MPAccessInterfaceType &mpAccess ()
    {
      assert( mpAccess_ );
      return *mpAccess_;
    }
    
  private:  
    // write data of DataImp& vector to object stream 
    // --writeBuffer 
    template< class DiscreteFunction >
    inline void writeBuffer ( const int link,
                              ObjectStreamType &str,
                              const DiscreteFunction &discreteFunction ) const
    {
      const IndexMapType &indexMap = sendIndexMap_[ linkRank_[link ] ];
      const int size = indexMap.size();

      typedef typename DiscreteFunction :: DofType DofType;
      enum { blockSize = DiscreteFunction :: 
              DiscreteFunctionSpaceType :: localBlockSize };
      // reserve buffer memory at once 
      str.reserve( str.size() + (size * blockSize * sizeof( DofType )) );

      // dirty hack to have faster access to stream 
      UnsaveObjectStream &os = (UnsaveObjectStream &)str;
      for( int i = 0; i < size; ++i )
      {
        // get dof block 
        typedef typename DiscreteFunction :: ConstDofBlockPtrType ConstDofBlockPtrType;
        ConstDofBlockPtrType blockPtr = discreteFunction.block( indexMap[i] );

        // write dof block to stream  
        for( int k = 0; k < blockSize; ++k )
        {
          os.writeUnsave( ((*blockPtr)[ k ]) );
        }
      }
    }

    // read data from object stream to DataImp& data vector 
    // --readBuffer 
    template< class DiscreteFunction, class Operation >
    inline void readBuffer ( const int link,
                             ObjectStreamType &str, 
                             DiscreteFunction &discreteFunction,
                             const Operation * ) const 
    {
      typedef typename DiscreteFunction :: DofType DofType;

      enum { blockSize = DiscreteFunction :: 
              DiscreteFunctionSpaceType :: localBlockSize };

      UnsaveObjectStream &os = (UnsaveObjectStream &)str;

      // get index map of rank belonging to link  
      const IndexMapType &indexMap = recvIndexMap_[ linkRank_[ link ] ];

      const int size = indexMap.size();
      for( int i = 0; i < size; ++i )
      {
        // get dof block 
        typedef typename DiscreteFunction :: DofBlockPtrType DofBlockPtrType;
        DofBlockPtrType blockPtr = discreteFunction.block( indexMap[i] );

        // read block 
        DofType value;
        for( int k = 0; k < blockSize; ++k )  
        {
          os.readUnsave( value );
          Operation :: apply( value, ((*blockPtr)[ k ]) );
        }
      }
    }
  };

  // --LinkBuilder 
  template< class Space >
  template< class LinkStorage, class IndexMapVector, InterfaceType CommInterface >
  class DependencyCache< Space > :: LinkBuilder
  : public CommDataHandleIF
    < LinkBuilder< LinkStorage, IndexMapVector, CommInterface >, int >
  {
  public:
    typedef LinkStorage LinkStorageType;

    typedef IndexMapVector IndexMapVectorType;

    typedef typename SpaceType :: BlockMapperType BlockMapperType; 

    typedef int DataType;

  protected:
    const int myRank_;
    const int mySize_; 
    
    LinkStorageType &linkStorage_; 

    IndexMapVectorType &sendIndexMap_;
    IndexMapVectorType &recvIndexMap_;

    const SpaceType &space_;
    const BlockMapperType &blockMapper_; 

  public:
    LinkBuilder ( LinkStorageType &linkStorage,
                  IndexMapVectorType &sendIdxMap,
                  IndexMapVectorType &recvIdxMap,
                  const SpaceType &space )
    : myRank_( space.grid().comm().rank() ),
      mySize_( space.grid().comm().size() ),
      linkStorage_( linkStorage ),
      sendIndexMap_( sendIdxMap ),
      recvIndexMap_( recvIdxMap ),
      space_( space ),
      blockMapper_( space.blockMapper() )
    {
    }

  protected:  
    void sendBackSendMaps() 
    {
      // create ALU communicator 
      MPAccessImplType mpAccess ( MPIHelper::getCommunicator() );

      // build linkage 
      mpAccess.removeLinkage();
      // insert new linkage 
      mpAccess.insertRequestSymetric( linkStorage_ );        
      // get destination ranks 
      std::vector<int> dest = mpAccess.dest();
      // get number of links 
      const int nlinks = mpAccess.nlinks();
      
      // create buffers 
      ObjectStreamVectorType osv( nlinks );

      //////////////////////////////////////////////////////////////
      //
      //  at this point complete send maps exsist on receiving side, 
      //  so send them back to sending side 
      //
      //////////////////////////////////////////////////////////////

      // write all send maps to buffer 
      for(int link=0; link<nlinks; ++link) 
      {
        sendIndexMap_[ dest[link] ].writeToBuffer( osv[link] );
      }

      // exchange data 
      osv = mpAccess.exchange( osv );

      // read all send maps from buffer 
      for(int link=0; link<nlinks; ++link) 
      {
        sendIndexMap_[ dest[link] ].readFromBuffer( osv[link] );
      }
    }
    
  public:  
    //! desctructor 
    ~LinkBuilder() 
    {
      sendBackSendMaps();
    }

    //! returns true if combination is contained 
    bool contains ( int dim, int codim ) const
    {
      return space_.contains( codim );
    }

    //! return whether we have a fixed size 
    bool fixedsize ( int dim, int codim ) const
    {
      return false;
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void gather ( MessageBuffer &buffer,
                  const Entity &entity ) const
    {
      // check whether we are a sending entity 
      const PartitionType myPartitionType = entity.partitionType();
      const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );

      // if we send data then send rank and dofs 
      if( send )
      {
        // send rank for linkage 
        buffer.write( myRank_ );

        // then send all dof numbers 
        const int numDofs = blockMapper_.numEntityDofs( entity );
        for( int i = 0; i < numDofs; ++i )
        {
          DataType idx = blockMapper_.mapEntityDofToGlobal( entity, i );
          buffer.write( idx );
        }
      }
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void scatter ( MessageBuffer &buffer,
                   const Entity &entity,
                   const size_t dataSize )
    {
      // if data size > 0 then other side is sender 
      if( dataSize > 0 ) 
      {
        // read rank of other side
        DataType rank;
        buffer.read( rank );  
        assert( (rank >= 0) && (rank < mySize_) );

        // check whether we are a sending entity 
        const PartitionType myPartitionType = entity.partitionType();
        const bool receive = EntityCommHelper< CommInterface > :: receive( myPartitionType );
        
        // insert rank of link into set of links
        linkStorage_.insert( rank );

        // read indices from stream 
        std::vector<int> indices( dataSize - 1 );
        for(size_t i=0; i<dataSize-1; ++i) 
        {
          buffer.read( indices[i] );  
        }

        // if we are a receiving entity 
        if( receive ) 
        {
          //////////////////////////////////////////////////////////
          //
          // Problem here: sending and receiving order might differ
          // Solution: sort all dofs after receiving order and send 
          // senders dofs back at the end 
          //
          //////////////////////////////////////////////////////////
          
          // if data has been send and we are receive entity 
          // then insert indices into send map of rank 
          sendIndexMap_[ rank ].insert( indices );
          
          // build local mapping for receiving of dofs 
          const int numDofs = blockMapper_.numEntityDofs( entity );
          indices.resize( numDofs );
          for( int i = 0; i < numDofs; ++i )
          {
            indices[ i ] = blockMapper_.mapEntityDofToGlobal( entity, i );
          }

          // insert receiving dofs 
          recvIndexMap_[ rank ].insert( indices );
        }
      }
    }

    //! return local dof size to be communicated 
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      const PartitionType myPartitionType = entity.partitionType();
      const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );
      return (send) ? (blockMapper_.numEntityDofs( entity ) + 1) : 0;
    }
  };



  //! object stream with unsave writing and reading
  template< class SpaceImp >
  class DependencyCache< SpaceImp > :: UnsaveObjectStream
  : public ALU3DSPACE ObjectStream 
  {
    typedef ALU3DSPACE ObjectStream BaseType;

  public:
    // create empty object stream 
    inline UnsaveObjectStream () : BaseType() {}
    // copy constructor taking object stream 
    inline UnsaveObjectStream (const ObjectStream & os) : BaseType(os) {}
    // copy constructor 
    inline UnsaveObjectStream (const UnsaveObjectStream & os) : BaseType(os) {}

    // write value to stream without testing size 
    template <class T> 
    inline void writeUnsave (const T & a)
    {
      T & val = *((T *) this->getBuff( this->_wb) );
      val = a;
      this->_wb += sizeof(T) ;
      assert( this->_wb <= this->_len );
      return ;
    } 
    
    // read value from stream without checking 
    template <class T>
    inline void readUnsave (T & a)
    {
      const T & val = *((const T *) this->getBuff(this->_rb) );
      a = val;
      this->_rb += sizeof(T);
      assert( this->_rb <= this->_wb ); 
      return ;
    }
  };



  template< class Space >
  inline void DependencyCache< Space > :: buildMaps ()
  {
    if( interface_ == InteriorBorder_All_Interface )
    {
      LinkBuilder< LinkStorageType, IndexMapVectorType,
                   InteriorBorder_All_Interface >
        handle( linkStorage_, sendIndexMap_, recvIndexMap_, space_ );
      buildMaps( handle );
    }
    else if( interface_ == InteriorBorder_InteriorBorder_Interface )
    {
      LinkBuilder< LinkStorageType, IndexMapVectorType,
                   InteriorBorder_InteriorBorder_Interface >
        handle( linkStorage_, sendIndexMap_, recvIndexMap_, space_ );
      buildMaps( handle );
    }
    else if( interface_ == All_All_Interface )
    {
      LinkBuilder< LinkStorageType, IndexMapVectorType, All_All_Interface >
        handle( linkStorage_, sendIndexMap_, recvIndexMap_, space_ );
      buildMaps( handle );
    }
    else
      DUNE_THROW( NotImplemented, "DependencyCache for the given interface has"
                                  " not been implemented, yet." );
#ifndef NDEBUG
    // checks that sizes of index maps are equal on sending and receiving proc 
    checkConsistency(); 
#endif
  }


  template< class Space >
  template< class LS, class IMV, InterfaceType CI >
  inline void DependencyCache< Space >
    :: buildMaps ( LinkBuilder< LS, IMV, CI > &handle )
  {
    linkStorage_.clear();
    for( int i = 0; i < mySize_; ++i )
    {
      recvIndexMap_[ i ].clear();
      sendIndexMap_[ i ].clear();
    }

    // make one all to all communication to build up communication pattern 
    gridPart_.communicate( handle, All_All_Interface , ForwardCommunication );

    // remove old linkage
    mpAccess().removeLinkage();
    // create new linkage
    mpAccess().insertRequestSymetric ( linkStorage_ );

    // get real rank numbers for each link 
    linkRank_ = mpAccess().dest();

    // remember number of links
    nLinks_ = mpAccess().nlinks();

    // resize buffer to number of links
    buffer_.resize( nLinks_ );

  }

  template< class Space >
  inline void DependencyCache< Space >
    :: checkConsistency ()  
  {
    /////////////////////////////
    // consistency check 
    /////////////////////////////

    // check that order and size are consistent 
    for(int l=0; l<nLinks_; ++l) 
    {
      buffer_[l].clear();
      const int sendSize = sendIndexMap_[ linkRank_[ l ] ].size();
      buffer_[l].write( sendSize );
      for(int i=0; i<sendSize; ++i) 
      {
        buffer_[l].write( i );
      }
    }

    // exchange data to other procs 
    buffer_ = mpAccess().exchange( buffer_  );
  
    // check that order and size are consistent 
    for(int l=0; l<nLinks_; ++l) 
    {
      const int recvSize = recvIndexMap_[ linkRank_[ l ] ].size();
      int sendedSize;
      buffer_[l].read( sendedSize );

      // compare sizes, must be the same 
      if( recvSize != sendedSize )
      {
        DUNE_THROW(InvalidStateException,"Sizes do not match!" << sendedSize << " o|r " << recvSize);
      }

      for(int i=0; i<recvSize; ++i) 
      {
        int idx;
        buffer_[l].read( idx );
        
        // ordering should be the same on both sides 
        if( i != idx ) 
        {
          DUNE_THROW(InvalidStateException,"Wrong ordering of send and recv maps!");
        }
      }
    }
  }

  template< class Space >
  template< class DiscreteFunction, class Operation >
  inline void DependencyCache< Space >
    :: exchange ( DiscreteFunction &discreteFunction,
                  const Operation *operation )
  {
    // on serial runs: do nothing 
    if( mySize_ <= 1 )
      return;

    // update cache 
    rebuild();
    
    const int numLinks = nlinks();

    // write buffers 
    for( int link = 0; link < numLinks; ++link )
    {
      buffer_[ link ].clear();
      writeBuffer( link, buffer_[ link ], discreteFunction );
    }

    // exchange data to other procs 
    buffer_ = mpAccess().exchange( buffer_  );
   
    // read buffers 
    for( int link = 0; link < numLinks; ++link )
      readBuffer( link, buffer_[ link ], discreteFunction, operation );
  }


  template< class Space >
  template< class DiscreteFunction >
  inline void DependencyCache< Space >
    :: writeBuffer ( ObjectStreamVectorType &osv,
                     const DiscreteFunction &discreteFunction ) const
  {
    const int numLinks = nlinks();
    // write buffers 
    for( int link = 0; link < numLinks; ++link )
      writeBuffer( link, osv[ link ], discreteFunction );
  }
 

  template< class Space >
  template< class DiscreteFunction, class Operation >
  inline void DependencyCache< Space >
    :: readBuffer ( ObjectStreamVectorType &osv,
                    DiscreteFunction &discreteFunction,
                    const Operation *operation ) const
  {
    const int numLinks = nlinks();
    // write buffers 
    for( int link = 0; link < numLinks; ++link )
      readBuffer( link, osv[ link ], discreteFunction, operation );
  }



  //! Key for CommManager singleton list 
  template <class SpaceImp>
  class CommManagerSingletonKey
  {
    const SpaceImp & space_;
    const InterfaceType interface_; 
    const CommunicationDirection dir_;
  public:
    //! constructor taking space 
    CommManagerSingletonKey(const SpaceImp & space,
                            const InterfaceType interface,
                            const CommunicationDirection dir) 
      : space_(space), interface_(interface), dir_(dir) {}

    //! copy constructor  
    CommManagerSingletonKey(const CommManagerSingletonKey & org) 
      : space_(org.space_), interface_(org.interface_), dir_(org.dir_) {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const CommManagerSingletonKey & otherKey) const
    {
      // mapper of space is singleton 
      return (&(space_.mapper()) == & (otherKey.space_.mapper()) );
    }

    //! return reference to index set 
    const SpaceImp & space() const { return space_; }
    //! return communication interface 
    const InterfaceType interface() const { return interface_; }
    //! return communication direction  
    const CommunicationDirection direction() const { return dir_; }
  };

  //! Factory class for SingletonList to tell how objects are created and
  //! how compared.
  template <class KeyImp, class ObjectImp>
  class CommManagerFactory
  {
    public:
    //! create new communiaction manager   
    static ObjectImp * createObject( const KeyImp & key )
    {
      return new ObjectImp(key.space(),
                           key.interface(),
                           key.direction());
    }
    
    //! delete comm manager  
    static void deleteObject( ObjectImp * obj )
    {
      delete obj; 
    }
  };

  //! Proxy class to DependencyCache which is singleton per space 
  template <class SpaceImp> 
  class CommunicationManager 
  {
    typedef CommunicationManager<SpaceImp> ThisType;
    
    // type of communication manager object which does communication 
    typedef DependencyCache<SpaceImp> DependencyCacheType;

    typedef CommManagerSingletonKey<SpaceImp> KeyType;
    typedef CommManagerFactory<KeyType, DependencyCacheType> FactoryType;

    typedef SingletonList< KeyType , DependencyCacheType , FactoryType > CommunicationProviderType;

    typedef SpaceImp SpaceType;
    const SpaceType & space_; 

    const KeyType key_;

    const int mySize_;

    typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType; 
    
    // is singleton per space 
    mutable DependencyCacheType &cache_;
    CommunicationManager(const ThisType& org);
  public:  
    //! constructor taking space and communication interface/direction 
    CommunicationManager(const SpaceType & space,
                         const InterfaceType interface,
                         const CommunicationDirection dir)
      : space_(space)
      , key_(space_,interface,dir)
      , mySize_(space_.grid().comm().size())
      , cache_(CommunicationProviderType::getObject(key_)) 
    {
    }

    //! constructor taking space and communication interface/direction 
    CommunicationManager(const SpaceType & space)
      : space_(space)
      , key_(space_,
             space.communicationInterface(),
             space.communicationDirection())
      , mySize_(space_.grid().comm().size())
      , cache_(CommunicationProviderType::getObject(key_)) 
    {
    }

    //! remove object comm
    ~CommunicationManager() 
    {
      CommunicationProviderType::removeObject(cache_);
    }

    //! return communication interface 
    InterfaceType communicationInterface() const {
      return cache_.communicationInterface();
    }

    //! return communcation direction 
    CommunicationDirection communicationDirection() const
    {
      return cache_.communicationDirection();
    }

    //! return time needed for last build 
    double buildTime() const { return cache_.buildTime(); }

    //! return time needed for last exchange  
    double exchangeTime() const { return cache_.exchangeTime(); }

    MPAccessInterfaceType& mpAccess() { return cache_.mpAccess(); }

    //! exchange discrete function to all procs we share data 
    //! using the copy operation 
    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) const
    {
      cache_.exchange( df, (DFCommunicationOperation :: Copy *) 0 );
    }
    
    //! exchange discrete function to all procs we share data 
    //! using the given operation 
    template <class DiscreteFunctionType, class OperationImp> 
    void exchange(DiscreteFunctionType & df, const OperationImp* ) const
    {
      cache_.exchange( df, (OperationImp*) 0 );
    }
    
    //! write given df to given buffer 
    template <class ObjectStreamVectorType, class DiscreteFunctionType> 
    void writeBuffer(ObjectStreamVectorType& osv, 
                     const DiscreteFunctionType & df) const 
    {
      cache_.writeBuffer(osv, df );
    }
    
    // read given df from given buffer 
    template <class ObjectStreamVectorType, class DiscreteFunctionType> 
    void readBuffer(ObjectStreamVectorType& osv, 
                    DiscreteFunctionType & df) const 
    {
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
        :: template CommDataHandle<DiscreteFunctionType> :: OperationType OperationType;
      cache_.readBuffer(osv, df , (OperationType *) 0 );
    }

    //! rebuild underlying cache if necessary 
    void rebuildCache () 
    {
      cache_.rebuild();
    }
  };

  //! Proxy class to DependencyCache which is singleton per space 
  class CommunicationManagerList  
  {
    //! communicated object interface 
    template <class MPAccessType, class ObjectStreamVectorType> 
    class DiscreteFunctionCommunicatorInterface  
    {
    protected:
      DiscreteFunctionCommunicatorInterface () {}
    public:
      virtual ~DiscreteFunctionCommunicatorInterface () {}

      virtual MPAccessType& mpAccess() = 0;
      virtual void writeBuffer(ObjectStreamVectorType&) const = 0;
      virtual void readBuffer(ObjectStreamVectorType&) = 0;
      virtual void rebuildCache () = 0;
    };
    
    //! communicated object implementation  
    template <class DiscreteFunctionImp,
              class MPAccessType,
              class ObjectStreamVectorType>
    class DiscreteFunctionCommunicator 
    : public DiscreteFunctionCommunicatorInterface<MPAccessType,ObjectStreamVectorType> 
    {
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
    
      // object to communicate 
      DiscreteFunctionType& df_;
      //! communicator manager 
      CommunicationManagerType comm_;
    public:  
      //! constructor taking disctete function 
      DiscreteFunctionCommunicator(DiscreteFunctionType& df) 
        : df_(df), comm_(df_.space())
      {}

      //! return ALUGrid communicator
      virtual MPAccessType& mpAccess() { return comm_.mpAccess(); }

      //! write discrete function to all buffers 
      virtual void writeBuffer(ObjectStreamVectorType& osv) const 
      {
        comm_.writeBuffer(osv,df_);
      }
      
      //! read discrete function from all buffers 
      virtual void readBuffer(ObjectStreamVectorType& osv)
      {
        comm_.readBuffer(osv,df_);
      }

      //! rebuild cache if grid changed 
      virtual void rebuildCache () 
      {
        comm_.rebuildCache ();
      }
    };

    // ALUGrid send/recv buffers 
    typedef ALU3DSPACE ObjectStream ObjectStreamType; 
    
    // type of buffer vector 
    typedef std::vector< ObjectStreamType > ObjectStreamVectorType;

    // type of ALUGrid Communicator 
    typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType; 

    // interface for communicated objects 
    typedef DiscreteFunctionCommunicatorInterface<MPAccessInterfaceType,ObjectStreamVectorType>
      CommObjInterfaceType;   
    
    // list of communicated objects 
    typedef std::list < CommObjInterfaceType * > CommObjListType; 
    CommObjListType objList_;

    // number of processors 
    int mySize_; 

    // copy constructor 
    CommunicationManagerList(const CommunicationManagerList&); 
  public:  
    //! constructor creating list of communicated objects 
    template <class CombinedObjectType>
    CommunicationManagerList(CombinedObjectType& cObj) : 
      mySize_( -1 )
    {
      // add all discrete functions containd in cObj to list 
      cObj.addToList(*this);
    }

    //! remove object comm
    ~CommunicationManagerList() 
    {
      // delete all entries 
      while( ! objList_.size() == 0 )
      {
        CommObjInterfaceType * obj = objList_.back();
        objList_.pop_back();
        delete obj;
      }
    }

    //! add one discrete function to the list 
    template <class DiscreteFunctionImp>
    void addToList(DiscreteFunctionImp &df)
    {
      // type of communication object 
      typedef DiscreteFunctionCommunicator<DiscreteFunctionImp,
                                           MPAccessInterfaceType,
                                           ObjectStreamVectorType> CommObj;
      CommObj * obj = new CommObj(df);
      objList_.push_back(obj);

      // if mySize wasn't set, set to number of processors
      if( mySize_ < 0 )
      {
        // get ALUGrid communicator 
        MPAccessInterfaceType& mpAccess = objList_.front()->mpAccess();

        // set number of processors 
        mySize_ = mpAccess.psize();
      } 
    }

    //! exchange the list of discrete functions between processes 
    //! only one communication is done here 
    void exchange() const
    {
      // if only one process, do nothing 
      if( mySize_ <= 1 ) return ;
        
      // exchange data 
      if(objList_.size() > 0)
      {
        typedef CommObjListType :: const_iterator iterator; 
        // rebuild cahce if grid has changed
        {
          iterator end = objList_.end();
          for(iterator it = objList_.begin(); it != end; ++it) 
          {
            (*it)->rebuildCache();
          }
        }
      
        // get ALUGrid communicator 
        MPAccessInterfaceType& mpAccess = objList_.front()->mpAccess();

        // create buffer 
        ObjectStreamVectorType osv( mpAccess.nlinks() );
        
        // fill buffers 
        iterator end  = objList_.end();
        for(iterator it = objList_.begin(); it != end; ++it) 
        {
          (*it)->writeBuffer(osv); 
        }
      
        // exchange data 
        osv = mpAccess.exchange (osv);

        // read buffers 
        for(iterator it = objList_.begin(); it != end; ++it) 
        {
          (*it)->readBuffer(osv); 
        }
      }
    }
  };
#endif 
  //@}
  
} // end namespace Dune 
#endif
