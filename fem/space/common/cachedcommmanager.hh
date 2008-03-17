#ifndef DUNE_CACHED_COMMUNICATION_MANAGER_HH
#define DUNE_CACHED_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/common/misc.hh>
#include <dune/common/mpihelper.hh>
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
   
    // for compatiblity with Robert's code
    enum { treatOverlapAsGhosts = true };
    
  protected:
    class CommunicationIndexMap;
    
    template< class LinkStorage, class IndexMapVector, InterfaceType CommInterface >
    class LinkBuilder;

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
    
  public:
    //! constructor taking space 
    DependencyCache ( const SpaceType &space,
                      const InterfaceType interface = InteriorBorder_All_Interface )
    : space_( space ),
      gridPart_( space_.gridPart() ),
      interface_( interface ),
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
      sequence_( -1 )
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
    // build linkage and index maps 
    inline void buildMaps ();

  protected:
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
      if( mySize_ <= 1 )
        return;

      // check whether grid has changed.
      if( sequence_ != space_.sequence() )
      {
        buildMaps();
        sequence_ = space_.sequence();
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
    template< class DiscreteFunction, class Operation >
    inline void readBuffer ( const int link,
                             ObjectStreamType &str, 
                             DiscreteFunction &discreteFunction,
                             const Operation *operation ) const 
    {
      typedef typename DiscreteFunction :: DofType DofType;

      enum { blockSize = DiscreteFunction :: 
              DiscreteFunctionSpaceType :: localBlockSize };

      UnsaveObjectStream &os = (UnsaveObjectStream &)str;
      
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



  // index map for send and receive data 
  template< class SpaceImp >
  class DependencyCache< SpaceImp > ::CommunicationIndexMap
  {
  protected:
    MutableArray< int > index_;

  public:
    //! constructor creating empty map
    CommunicationIndexMap() : index_(0) 
    {
      index_.setMemoryFactor( 1.1 );
    }

  private:
    // prohibit copying
    CommunicationIndexMap( const CommunicationIndexMap & );

  public:
    //! reserve memory 
    void reserve( int size ) 
    {
      // resize array, memory factor will be used 
      index_.resize( size );
    }

    //! clear index map 
    void clear() 
    {
      // resize 0 will free memory 
      index_.resize( 0 );
    }

    //! append index vector with idx 
    void insert( const std::vector<int> & idx )
    {
      const int size = idx.size();
      int count = index_.size();
      // reserve memory 
      reserve( count + size );
      assert( index_.size() == (count+size));
      // copy indices to index vector 
      for(int i=0; i<size; ++i, ++count) 
      { 
        assert( idx[i] >= 0 );
        index_[count] = idx[i]; 
      }
    }

    //! return index map for entry i
    const int operator [] (int i) const 
    {
      assert( i >= 0 );
      assert( i < (int) index_.size());
      return index_[i];
    }

    //! return size of map
    int size () const { return index_.size(); }

    //! print  map for debugging only 
    void print(std::ostream & s, int rank) const 
    {
      const int size = index_.size();
      s << "Start print: size = " << size << std::endl;
      for(int i=0; i<size; ++i) 
      {
        s<< rank << " idx["<<i<<"] = " << index_[i] << std::endl;
      }
      s << "End of Array" << std :: endl;
    }
  };



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
    {}

    bool contains ( int dim, int codim ) const
    {
      return space_.contains( codim );
    }

    bool fixedsize ( int dim, int codim ) const
    {
      return true;
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void gather ( MessageBuffer &buffer,
                  const Entity &entity ) const
    {
      // send rank 
      buffer.write( myRank_ );
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void scatter ( MessageBuffer &buffer,
                   const Entity &entity,
                   const size_t dataSize )
    {
      // build local mapping 
      const int numDofs = blockMapper_.numEntityDofs( entity );
      std :: vector< int > indices( numDofs );
      for( int i = 0; i < numDofs; ++i )
      {
        indices[ i ] = blockMapper_.mapEntityDofToGlobal( entity, i );
      }

      const PartitionType p = entity.partitionType();
      const bool send = EntityCommHelper< CommInterface > :: send( p );
      const bool receive = EntityCommHelper< CommInterface > :: receive( p );
      
      // read links and insert to mappings 
      for( size_t i = 0; i < dataSize; ++i )
      {
        // read rank of other side
        DataType value;
        buffer.read( value );  
        assert( (value >= 0) && (value < mySize_) );
        
        // insert rank of link into set of links
        linkStorage_.insert( value );

        // if we are on a send entity, insert into send cache 
        if( send )
          sendIndexMap_[ value ].insert( indices );

        // if we are on a receive entity, insert into receive cache 
        if( receive ) 
          recvIndexMap_[ value ].insert( indices );
      }
    }

    //! return local dof size to be communicated 
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      return 1; 
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

    // do communication to build up linkage
    if( treatOverlapAsGhosts && (gridPart_.grid().overlapSize( 0 ) > 0) )
    {
      // case of YaspGrid, where Overlap is treated as ghost 
      gridPart_.communicate
        ( handle, InteriorBorder_All_Interface, ForwardCommunication );
      gridPart_.communicate
        ( handle, InteriorBorder_All_Interface, BackwardCommunication );
    }
    else 
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
  public:
    //! constructor taking space 
    CommManagerSingletonKey(const SpaceImp & space) : space_(space) {}
    //! copy constructor  
    CommManagerSingletonKey(const CommManagerSingletonKey & org) : space_(org.space_) {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const CommManagerSingletonKey & otherKey) const
    {
      // mapper of space is singleton 
      return (&(space_.mapper()) == & (otherKey.space_.mapper()) );
    }

    //! return reference to index set 
    const SpaceImp & space() const { return space_; }
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
      return new ObjectImp(key.space());
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
    DependencyCacheType & cache_;
    CommunicationManager(const ThisType& org);
  public:  
    //! constructor taking space 
    CommunicationManager(const SpaceType & space) 
      : space_(space)
      , key_(space_)
      , mySize_(space_.grid().comm().size())
      , cache_(CommunicationProviderType::getObject(key_)) 
    {
    }

    //! remove object comm
    ~CommunicationManager() 
    {
      CommunicationProviderType::removeObject(cache_);
    }

    MPAccessInterfaceType& mpAccess() { return cache_.mpAccess(); }

    //! exchange discrete function to all procs we share data 
    //! using the copy operation 
    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) 
    {
      cache_.exchange( df, (DFCommunicationOperation :: Copy *) 0 );
    }
    
    //! exchange discrete function to all procs we share data 
    //! using the given operation 
    template <class DiscreteFunctionType, class OperationImp> 
    void exchange(DiscreteFunctionType & df, const OperationImp* ) 
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
    void exchange() 
    {
      // if only one process, do nothing 
      if( mySize_ <= 1 ) return ;
        
      // exchange data 
      if(objList_.size() > 0)
      {
        typedef CommObjListType :: iterator iterator; 
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
