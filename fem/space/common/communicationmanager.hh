#ifndef DUNE_COMMUNICATION_MANAGER_HH
#define DUNE_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/datahandleif.hh>
#if HAVE_ALUGRID
// inlcude alugrid to have to communicator class from ALUGrid 
#include <dune/grid/alugrid.hh>
#endif

//- Dune-fem includes 
#include <dune/fem/function/common/dfcommunication.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune { 
  
/** @addtogroup Communication Communication 
    @{
**/
  
  //! \brief Default CommunicationManager class just using the grids communicate
  //! method 
  template <class SpaceImp, 
            class OperationImp = DFCommunicationOperation::Copy >
  class DefaultCommunicationManager 
  {
    typedef SpaceImp SpaceType; 
    typedef typename SpaceType :: GridPartType GridPartType; 

    // gridPart for communication 
    const GridPartType & gridPart_; 

    DefaultCommunicationManager(const DefaultCommunicationManager &);
  public:
    //! constructor taking space, but here only storing gridPart for
    //! communication
    DefaultCommunicationManager(const SpaceType & space)
      : gridPart_(space.gridPart()) 
    {}

    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) 
    {
      // if serial run, just return   
      if(gridPart_.grid().comm().size() <= 1) return;
     
      // get data handler 
      typedef DiscreteFunctionCommunicationHandler<DiscreteFunctionType,OperationImp>
           DataHandleType;
      DataHandleType dataHandle(df);

      // communicate data 
      gridPart_.communicate( dataHandle, InteriorBorder_All_Interface , ForwardCommunication);
    }
  };
  
  // only if ALUGrid found and was build for parallel runs 
#if HAVE_ALUGRID && ALU3DGRID_PARALLEL 
  //! class to build up a map of all dofs of entities to be exchanged
  //! during a communication procedure. This speeds up the communication
  //! procedure, because no grid traversal is necessary to exchange data.
  //! this class is singleton for different discrete function space,
  //! because the dof mapping is always the same.
  template <class SpaceImp> 
  class DependencyCache 
  {
    // index map for send and receive data 
    class CommunicationIndexMap
    {
      //MutableArray<int> index_;  
      std::vector<int> index_;
      CommunicationIndexMap(const CommunicationIndexMap&);
    public:
      //! constructor creating empty map
      CommunicationIndexMap() : index_(0) {}

      //! reserve memory 
      void reserve( int size ) 
      {
        index_.reserve( size );
      }

      //! clear index map 
      void clear() 
      {
        index_.clear();
        index_.resize( 0 );
      }

      //! append index vector with idx 
      void insert( const std::vector<int> & idx )
      {
        const int size = idx.size();
        const int newSize = index_.size() + idx.size();
        reserve( newSize );
        for(int i=0; i<size; ++i) 
        { 
          assert( idx[i] >= 0 );
          index_.push_back( idx[i] ); 
        }
        assert( (int) index_. size () == newSize );
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
        s << "End of Array\n";
      }
    };

    template <class LinkStorageImp, class IndexMapVectorType, 
              class SpaceType> 
    class LinkBuilder
     : public CommDataHandleIF<
       LinkBuilder< LinkStorageImp , IndexMapVectorType , SpaceType > ,
       int >
    {
    public:
      typedef int DataType;

      const int myRank_;
      const int mySize_; 
      
      typedef LinkStorageImp LinkStorageType;
      // discrete function to communicate 
      mutable LinkStorageType & linkStorage_; 

      IndexMapVectorType & sendIndexMap_;
      IndexMapVectorType & recvIndexMap_;

      const SpaceType & space_;

    public:
      LinkBuilder(LinkStorageType & linkStorage,
            IndexMapVectorType & sendIdxMap, 
            IndexMapVectorType & recvIdxMap,
            const SpaceType & space)
        : myRank_(space.grid().comm().rank()) 
        , mySize_(space.grid().comm().size())
        , linkStorage_(linkStorage)
        , sendIndexMap_(sendIdxMap)
        , recvIndexMap_(recvIdxMap)
        , space_(space)
      {
      }

      bool contains (int dim, int codim) const
      {
        return (codim == 0);
      }

      bool fixedsize (int dim, int codim) const
      {
        return true;
      }

      //! read buffer and apply operation 
      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& en) const
      {
        // send rank 
        buff.write( myRank_ );
      }

      //! read buffer and apply operation 
      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
      {
        // assert that we get the same size of data we sent 
        assert( n == 1 );
        
        // build local mapping 
        const int numDofs = space_.baseFunctionSet(en).numBaseFunctions(); 
        std::vector<int> indices(numDofs);
        // copy numDofs 
        for(int i=0; i<numDofs; ++i) 
        {
          indices[i] = space_.mapToGlobal( en , i ); 
        }

        // Interior entities belong to send area and other entities, i.e.
        // Overlap and Ghost, belong to receive area 
        const bool interiorEn = (en.partitionType() == InteriorEntity); 
        
        // read links and insert mapping 
        DataType val;
        
        // read rank of other side 
        buff.read( val );  
        
        // check that value of rank is within valid range 
        assert( val < mySize_ );
        assert( val >= 0 );
      
        // insert rank of link into set of links
        linkStorage_.insert( val );

        if(interiorEn) 
        {
          sendIndexMap_[val].insert( indices );
        }
        else 
        {
          recvIndexMap_[val].insert( indices );
        }
      }

      //! return local dof size to be communicated 
      template<class EntityType>
      size_t size (const EntityType& en) const
      {
        // size of data 
        return 1; 
      }
    };

    //! type of discrete function space 
    typedef SpaceImp SpaceType; 
    //! type of grid part 
    typedef typename SpaceType :: GridPartType GridPartType; 

    const SpaceType & space_; 
    const GridPartType & gridPart_; 

    const int myRank_;
    const int mySize_; 
    
    typedef std::set<int>  LinkStorageType;
    // type of set of links 
    LinkStorageType linkStorage_; 

    // type of communication indices 
    typedef CommunicationIndexMap IndexMapType;
    // type of IndexMapVector 
    typedef IndexMapType* IndexMapVectorType;
    IndexMapType* recvIndexMap_;
    IndexMapType* sendIndexMap_;

    // vector containing the links of this process 
    std::vector < int > linkRank_;

    // ALUGrid send/recv buffers 
    typedef ALU3DSPACE ObjectStream ObjectStreamType; 
    
    // type of communicator 
    typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType; 
    // type of communication implementation 
    typedef ALU3DSPACE MpAccessMPI   MPAccessImplType; 

    // ALUGrid communicatior Class 
    MPAccessInterfaceType * mpAccess_;

    //! type of communication buffer vector 
    typedef std::vector< ObjectStreamType > ObjectStreamVectorType;
    //! communication buffers 
    ObjectStreamVectorType buffer_;

    //! number of links 
    int nLinks_;

    //! know grid sequence number 
    int sequence_; 
    
    // do not copy this class 
    DependencyCache(const DependencyCache &);
  public:
    //! constructor taking space 
    DependencyCache(const SpaceType & space)
      : space_(space) , gridPart_(space_.gridPart()) 
      , myRank_(gridPart_.grid().comm().rank())
      , mySize_(gridPart_.grid().comm().size())
      , linkStorage_()
      , recvIndexMap_(new IndexMapType[mySize_])
      , sendIndexMap_(new IndexMapType[mySize_])
      , linkRank_()
      // create mpAccess with world communicator  
      // only when size > 1 
      , mpAccess_( (mySize_ > 1) ? 
          (new MPAccessImplType( MPIHelper::getCommunicator() )) : 0)
      , nLinks_(0)
      , sequence_(-1)
    {
    }

    //! destrcutor removeing mpAccess 
    ~DependencyCache()
    {
      delete mpAccess_; mpAccess_ = 0;
      delete [] sendIndexMap_; sendIndexMap_ = 0;
      delete [] recvIndexMap_; recvIndexMap_ = 0;
    }

  public:
    // build linkage and index maps 
    void buildMaps() 
    {
      linkStorage_.clear();
      for(int i=0; i<mySize_; ++i)
      {
        recvIndexMap_[i].clear();
        sendIndexMap_[i].clear();
      }

      // type of data handler 
      typedef LinkBuilder< LinkStorageType , IndexMapVectorType , SpaceType > LinkBuilderHandleType; 
      LinkBuilderHandleType handle( linkStorage_, sendIndexMap_, recvIndexMap_ , space_ );

      // do communication to build up linkage 
      if( gridPart_.grid().overlapSize(0) > 0 ) 
      {
        // case of YaspGrid, where Overlap is treated as ghost 
        gridPart_.communicate( handle, InteriorBorder_All_Interface , ForwardCommunication );
        gridPart_.communicate( handle, InteriorBorder_All_Interface , BackwardCommunication );
      }
      else 
      {
        // case of ALUGrid, where we have only the interior ghost situation 
        gridPart_.communicate( handle, All_All_Interface , ForwardCommunication );
      }

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

    //! return number of links 
    int nlinks () const { return nLinks_; }

    //! check if grid has changed and rebuild cache if necessary 
    void rebuild() 
    {
      if(sequence_ != space_.sequence()) 
      {
        buildMaps();

        // store actual sequence number 
        sequence_ = space_.sequence();
      }
    }
      
    //! exchange data of discrete function 
    template<class DiscreteFunctionType, class OperationImp>
    void exchange(DiscreteFunctionType& df, const OperationImp *)
    {
      typedef typename DiscreteFunctionType :: DofType DofType;
      // if serial run, just return   
      if(mySize_ <= 1) return;
       
      // check if rebuild is needed and update cache 
      rebuild();
      
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        // reset buffers, keeps memory  
        buffer_[l].clear();

        writeBuffer( l , buffer_[l] , df.leakPointer(), (DofType *) 0 );
      }

      // exchange data to other procs 
      buffer_ = mpAccess().exchange( buffer_ );
     
      // read buffers 
      for(int l=0; l<links; ++l) 
      {
        readBuffer( l , buffer_[l] , df.leakPointer() ,
                    (DofType *) 0, (OperationImp *) 0 );
      }
    }

    //! write data of discrete function to buffer 
    template<class DiscreteFunctionType>
    void writeBuffer(ObjectStreamVectorType & osv,
                     DiscreteFunctionType& df) const 
    {
      typedef typename DiscreteFunctionType :: DofType DofType;
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        writeBuffer( l , osv[l] , df.leakPointer(), (DofType *) 0);
      }
    }
    
    //! read data of discrete function from buffer  
    template<class DiscreteFunctionType, class OperationImp>
    void readBuffer(ObjectStreamVectorType & osv,
                    DiscreteFunctionType& df, const OperationImp *o) const 
    {
      typedef typename DiscreteFunctionType :: DofType DofType;
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        readBuffer( l, osv[l] , df.leakPointer(), 
                    (DofType*) 0, (OperationImp *) 0); 
      }
    }
    
    //! return reference to mpAccess object
    MPAccessInterfaceType& mpAccess() 
    {   
      assert( mpAccess_ );
      return *mpAccess_;
    }
    
  private:  
    // write data of DataImp& vector to object stream 
    template <class DataImp, class DofType> 
    void writeBuffer(const int link, 
                     ObjectStreamType & os, 
                     const DataImp& data,
                     const DofType* ) const 
    {
      const IndexMapType& indexMap = sendIndexMap_[ linkRank_ [link ] ]; 
      const int size = indexMap.size();

      // reserve buffer memory at once 
      os.reserve( size * sizeof(DofType) );

      DofType val = 0.0;
      for(int i=0; i<size; ++i)
      {
        val = data[ indexMap[i] ];
        os.write( val );
      }
    }

    // read data from object stream to DataImp& data vector 
    template <class DataImp, class DofType, class OperationImp> 
    void readBuffer(const int link, 
                    ObjectStreamType & os, 
                    DataImp& data, const DofType*,
                    const OperationImp *) const 
    {
      const IndexMapType& indexMap = recvIndexMap_[ linkRank_ [link ] ]; 
      DofType val;
      const int size = indexMap.size();
      for(int i=0; i<size; ++i)
      {
        os.read( val );
        // apply operation 
        OperationImp::apply(val , data[ indexMap[i] ] );
      }
    }
    
    // write data of double* vector to object stream 
    void writeBuffer(const int link, 
                     ObjectStreamType & os, 
                     const double* data,
                     const double* ) const 
    {
      const IndexMapType& indexMap = sendIndexMap_[ linkRank_ [link ] ]; 
      const int size = indexMap.size();

      // reserve buffer memory at once 
      os.reserve( size * sizeof(double) );

      double val = 0.0;
      for(int i=0; i<size; ++i)
      {
        val = data[ indexMap[i] ];
        os.write( val );
      }
    }

    // read data from object stream to double* data vector 
    template <class OperationImp> 
    void readBuffer(const int link, 
                    ObjectStreamType & os, 
                    double* data, const double*,
                    const OperationImp *) const 
    {
      const IndexMapType& indexMap = recvIndexMap_[ linkRank_ [link ] ]; 
      double val;
      const int size = indexMap.size();
      for(int i=0; i<size; ++i)
      {
        os.read( val );
        // apply operation 
        OperationImp::apply(val , data[ indexMap[i] ] );
      }
    }
    
  };

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
  template <class SpaceImp, 
            class OperationImp = DFCommunicationOperation::Copy >
  class CommunicationManager 
  {
    typedef CommunicationManager<SpaceImp,OperationImp> ThisType;
    
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

    //! exchange discrete function to all procs we share data with 
    //! by using given OperationImp when receiving data from other procs 
    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) 
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
      cache_.readBuffer(osv, df , (OperationImp *) 0 );
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
              class ObjectStreamVectorType,
              class OperationImp>
    class DiscreteFunctionCommunicator 
    : public DiscreteFunctionCommunicatorInterface<MPAccessType,ObjectStreamVectorType> 
    {
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef CommunicationManager<DiscreteFunctionSpaceType,OperationImp> CommunicationManagerType; 
    
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

    // copy constructor 
    CommunicationManagerList(const CommunicationManagerList&); 
  public:  
    //! constructor creating list of communicated objects 
    template <class CombinedObjectType>
    CommunicationManagerList(CombinedObjectType& cObj) 
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
                                           ObjectStreamVectorType,
                                           DFCommunicationOperation::Copy> CommObj;
      CommObj * obj = new CommObj(df);
      objList_.push_back(obj);
    }

    //! exchange the list of discrete functions between processes 
    //! only one communication is done here 
    void exchange() 
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
      
      // exchange data 
      if(objList_.size() > 0)
      {
        // get ALUGrid communicator 
        MPAccessInterfaceType& mpAccess = objList_.front()->mpAccess();

        // if only one process, do nothing 
        if( mpAccess.psize() <= 1 ) return ;
        
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

#else 
  // if no ALUGrid found, supply default implementation 
#ifndef NDEBUG 
#if HAVE_MPI == 0
  #warning "HAVE_MPI == 0, therefore default CommunicationManager is used!"
#elif !ALU3DGRID_PARALLEL 
  #warning "No Parallel ALUGrid found, using default CommunicationManager!"
#endif 
#endif
  //! \brief use Default CommunicationManager as Communication Manager 
  template <class SpaceImp, 
            class OperationImp = DFCommunicationOperation::Copy >
  class CommunicationManager 
  : public DefaultCommunicationManager<SpaceImp,OperationImp> 
  {
    typedef DefaultCommunicationManager<SpaceImp,OperationImp> BaseType;
    CommunicationManager(const CommunicationManager &);
  public:
    //! constructor taking space, but here only storing gridPart for
    //! communication
    CommunicationManager(const SpaceImp & space) 
      : BaseType(space) 
    {}
  };

  //! Proxy class to DependencyCache which is singleton per space 
  class CommunicationManagerList  
  {
    //! communicated object interface 
    class DiscreteFunctionCommunicatorInterface  
    {
    protected:
      DiscreteFunctionCommunicatorInterface () {}
    public:
      virtual ~DiscreteFunctionCommunicatorInterface () {}
      virtual void exchange () = 0;
    };
    
    //! communicated object implementation  
    template <class DiscreteFunctionImp>
    class DiscreteFunctionCommunicator 
    : public DiscreteFunctionCommunicatorInterface 
    {
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      
      // operation to perform 
      typedef typename DFCommunicationOperation :: Copy OperationType;
      
      typedef CommunicationManager<DiscreteFunctionSpaceType,OperationType> CommunicationManagerType; 
    
      typedef DiscreteFunctionCommunicationHandler<DiscreteFunctionType,OperationType> DataHandleType;
     
      DiscreteFunctionType& df_;
      CommunicationManagerType comm_;
    public:  
      //! constructor taking disctete function 
      DiscreteFunctionCommunicator(DiscreteFunctionType& df) 
        : df_(df), comm_(df_.space())
      {
      }

      // exchange discrete function 
      void exchange () 
      {
        comm_.exchange(df_);
      }
    };

    typedef DiscreteFunctionCommunicatorInterface CommObjIFType;
    typedef std::list < DiscreteFunctionCommunicatorInterface * > CommObjListType;
    CommObjListType objList_;

    CommunicationManagerList(const CommunicationManagerList&); 
  public: 
    //! constructor
    template <class CombinedObjectType>
    CommunicationManagerList(CombinedObjectType& cObj) 
    {
      cObj.addToList(*this);
    }

    //! remove object comm
    ~CommunicationManagerList() 
    {
      // delete all entries 
      while( ! objList_.size() == 0 )
      {
        CommObjIFType * obj = objList_.back();
        objList_.pop_back();
        delete obj;
      }
    }

    //! add discrete function to communication list 
    template <class DiscreteFunctionImp> 
    void addToList(DiscreteFunctionImp &df)
    {
      typedef DiscreteFunctionCommunicator<DiscreteFunctionImp> CommObjType;
      CommObjType* obj = new CommObjType(df);
      objList_.push_back(obj);
    }

    //! exchange discrete function to all procs we share data with 
    //! by using given OperationImp when receiving data from other procs 
    void exchange() 
    {
      typedef CommObjListType :: iterator iterator; 
      {
        iterator end = objList_.end();
        for(iterator it = objList_.begin(); it != end; ++it) 
        {
          (*it)->exchange();
        }
      }
    }
  };

  // end toggle AULGrid yes/no
#endif 
  //@}
  
} // end namespace Dune 
#endif
