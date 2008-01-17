#ifndef DUNE_CACHED_COMMUNICATION_MANAGER_HH
#define DUNE_CACHED_COMMUNICATION_MANAGER_HH

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
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune { 
  
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
  template <class SpaceImp> 
  class DependencyCache 
  {
    // index map for send and receive data 
    class CommunicationIndexMap
    {
      MutableArray<int> index_;
      CommunicationIndexMap(const CommunicationIndexMap&);
    public:
      //! constructor creating empty map
      CommunicationIndexMap() : index_(0) 
      {
        index_.setMemoryFactor(1.1);
      }

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
      template <int dummy, int codim> 
      struct CheckInterior 
      {
        inline static bool check(const PartitionType p) 
        {
          DUNE_THROW(NotImplemented,"Method not implemented!");
          return true; 
        }
      };

      // codim 0 specialization 
      template <int dummy> 
      struct CheckInterior<dummy,0>
      {
        inline static bool check(const PartitionType p) 
        {
          return (p == InteriorEntity); 
        }
      };
      
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
      typedef typename SpaceType :: MapperType MapperType; 
      const MapperType& mapper_; 

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
        , mapper_(space.mapper())
      {
      }

      bool contains (int dim, int codim) const
      {
        return space_.contains(codim);
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
      void scatter (MessageBufferImp& buff, 
                    const EntityType& en, 
                    const size_t dataSize)
      {
        // build local mapping 
        const int numDofs = mapper_.numEntityDofs( en );
        std::vector<int> indices(numDofs);

        // copy numDofs 
        for(int i=0; i<numDofs; ++i) 
        {
          indices[i] = mapper_.mapEntityDofToGlobal( en , i ); 
        }

        // Interior entities belong to send area and other entities, i.e.
        // Overlap and Ghost, belong to receive area 
        const bool interiorEn = 
          CheckInterior<-1,EntityType :: codimension>::check(en.partitionType());
        
        // read links and insert to mappings 
        for(size_t i=0; i<dataSize; ++i) 
        {
          // create data type 
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
      }

      //! return local dof size to be communicated 
      template<class EntityType>
      size_t size (const EntityType& en) const
      {
        // size of data 
        return 1; 
      }
    };

    //! object stream with unsave writing and reading 
    class UnsaveObjectStream : public ALU3DSPACE ObjectStream 
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
    typedef ALU3DSPACE ObjectStream   ObjectStreamType; 
    
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
      //, mpAccess_( (mySize_ > 1) ? 
      //    (new MPAccessImplType( MPIHelper::getCommunicator() )) : 0)
      , mpAccess_( new MPAccessImplType( MPIHelper::getCommunicator() ) )
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
      // only in parallel we have to do something 
      if( mySize_ <= 1 ) return ;

      // check whether grid has changed. 
      if(sequence_ != space_.sequence()) 
      {
        buildMaps();

        // store actual sequence number 
        sequence_ = space_.sequence();
      }
    }
      
    //! exchange data of discrete function 
    template<class DiscreteFunctionType, class OperationImp >
    void exchange(DiscreteFunctionType& discreteFunction, 
                  const OperationImp *)
    {
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

        writeBuffer( l , buffer_[l] , discreteFunction );
      }

      // exchange data to other procs 
      buffer_ = mpAccess().exchange( buffer_ );
     
      // read buffers 
      for(int l=0; l<links; ++l) 
      {
        readBuffer( l , buffer_[l] , discreteFunction, (OperationImp *) 0 );
      }
    }

    //! write data of discrete function to buffer 
    template<class DiscreteFunctionType>
    void writeBuffer(ObjectStreamVectorType & osv,
                     const DiscreteFunctionType& discreteFunction) const 
    {
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        writeBuffer( l , osv[l] , discreteFunction);
      }
    }
    
    //! read data of discrete function from buffer  
    template<class DiscreteFunctionType, class OperationImp>
    void readBuffer(ObjectStreamVectorType & osv,
                    DiscreteFunctionType& discreteFunction, 
                    const OperationImp *o) const 
    {
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        readBuffer( l, osv[l] , discreteFunction, (OperationImp *) 0); 
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
    template <class DiscreteFunctionImp> 
    void writeBuffer(const int link, 
                     ObjectStreamType & str, 
                     const DiscreteFunctionImp& discreteFunction) const 
    {
      const IndexMapType& indexMap = sendIndexMap_[ linkRank_ [link ] ]; 
      const int size = indexMap.size();

      typedef typename DiscreteFunctionImp :: DofType DofType;
      // reserve buffer memory at once 
      str.reserve( str.size() + (size * sizeof(DofType)) );

      // dirty hack to have faster access to stream 
      UnsaveObjectStream& os = (UnsaveObjectStream &) str;
      for(int i=0; i<size; ++i)
      {
        os.writeUnsave( discreteFunction.dof( indexMap[i] ) );
      }
    }

    // read data from object stream to DataImp& data vector 
    template <class DiscreteFunctionImp, class OperationImp> 
    void readBuffer(const int link, 
                    ObjectStreamType & str, 
                    DiscreteFunctionImp& discreteFunction,
                    const OperationImp *) const 
    {
      UnsaveObjectStream& os = (UnsaveObjectStream &) str;
      const IndexMapType& indexMap = recvIndexMap_[ linkRank_ [link ] ]; 
      typedef typename DiscreteFunctionImp :: DofType DofType;
      DofType val;
      const int size = indexMap.size();
      for(int i=0; i<size; ++i)
      {
        os.readUnsave( val );
        // apply operation 
        OperationImp::apply(val , discreteFunction.dof( indexMap[i] ) );
      }
    }
    
    // write data of double* vector to object stream 
    template <class T>
    void writeBuffer(const int link, 
                     ObjectStreamType & str, 
                     const T* data,
                     const T* ) const 
    {
      const IndexMapType& indexMap = sendIndexMap_[ linkRank_ [link ] ]; 
      const int size = indexMap.size();

      // reserve buffer memory at once 
      str.reserve( str.size() + (size * sizeof(T)) );

      // dirty hack to have faster access to stream 
      UnsaveObjectStream& os = (UnsaveObjectStream &) str;
      T val = 0;
      for(int i=0; i<size; ++i)
      {
        val = data[ indexMap[i] ];
        os.writeUnsave( val );
      }
    }

    // read data from object stream to double* data vector 
    template <class OperationImp, class T> 
    void readBuffer(const int link, 
                    ObjectStreamType & str, 
                    T* data, const T*,
                    const OperationImp *) const 
    {
      const IndexMapType& indexMap = recvIndexMap_[ linkRank_ [link ] ]; 

      // dirty hack to have faster access to stream 
      UnsaveObjectStream& os = (UnsaveObjectStream &) str;

      T val;
      const int size = indexMap.size();
      for(int i=0; i<size; ++i)
      {
        os.readUnsave( val );
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
