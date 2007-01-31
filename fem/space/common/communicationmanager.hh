#ifndef DUNE_COMMUNICATION_MANAGER_HH
#define DUNE_COMMUNICATION_MANAGER_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/grid/common/datahandleif.hh>
#if HAVE_ALUGRID
// inlcude alugrid to have to communicator class from ALUGrid 
#include <dune/grid/alugrid.hh>
#endif

//- Dune-fem includes 
#include <dune/fem/function/common/dfcommunication.hh>
#include <dune/fem/space/common/singletonlist.hh>

namespace Dune { 

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
      std::vector<int> index_;  
    public:
      CommunicationIndexMap() {}

      void reserve( int size ) 
      {
        index_.reserve( size );
      }

      void clear() 
      {
        index_.resize( 0 );
      }

      void insert( const std::vector<int> & idx )
      {
        const int size = idx.size();
        const int newSize = index_.size() + size;
        reserve( newSize );
        for(int i=0; i<size; ++i) 
        { 
          index_.push_back( idx[i] ); 
        }
        assert( (int) index_. size () == newSize );
      }

      const int & operator [] (int i) const 
      {
        assert( i >= 0 );
        assert( i < (int) index_.size());
        return index_[i];
      }

      int size () const { return index_.size(); }

      void print(std::ostream & s, int rank) const 
      {
        const int size = index_.size();
        s << "Start print: size = " << size << std::endl;
        for(int i=0; i<size; ++i) 
        {
          s<< rank << " idx = " << index_[i] << std::endl;
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
        for(int i=0; i<numDofs; ++i) 
        {
          indices[i] = space_.mapToGlobal( en , i ); 
        }
      
        // Interior entities belong to send area and other entities, i.e.
        // Overlap and Ghost, belong to receive area 
        const bool interiorEn = (en.partitionType() == InteriorEntity); 
        
        // read links and insert mapping 
        DataType val;
        
        buff.read( val );  
        // size of sendIndexMap_ is equal to number of procs 
        assert( val <= (int) sendIndexMap_.size() );
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
    LinkStorageType linkStorage_; 

    typedef CommunicationIndexMap IndexMapType;
    typedef std::vector<IndexMapType> IndexMapVectorType;
    // index maps for receiving data 
    IndexMapVectorType recvIndexMap_;
    // index maps for sending data 
    IndexMapVectorType sendIndexMap_;

    std::vector < int > linkRank_;

    // ALUGrid send/recv buffers 
    typedef ALU3DSPACE ObjectStream ObjectStreamType; 
    // ALUGrid communicatior Class 
    ALU3DSPACE MpAccessMPI mpAccess_;

    //! communication buffers 
    std::vector< ObjectStreamType > buffer_;

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
      , recvIndexMap_(mySize_)
      , sendIndexMap_(mySize_)
      , linkRank_()
      // create mpAccess with communicator from grid 
      , mpAccess_(gridPart_.grid().comm().operator MPI_Comm())
      , nLinks_(0)
      , sequence_(-1)
    {
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
      mpAccess_.removeLinkage(); 
      // create new linkage 
      mpAccess_.insertRequestSymetric ( linkStorage_ );

      // get real rank numbers for each link 
      linkRank_ = mpAccess_.dest();

      // remember number of links 
      nLinks_ = mpAccess_.nlinks();

      buffer_.resize( nLinks_ );
      //std::cout << "Build dependency cache. Size = " << sendIndexMap_[linkRank_[0 ]].size() <<"\n";
    }

    //! return number of links 
    int nlinks () const { return nLinks_; }

    //! exchange data of discrete function 
    template<class DiscreteFunctionType, class OperationImp>
    void exchange(DiscreteFunctionType& df, const OperationImp *)
    {
      // if serial run, just return   
      if(mySize_ <= 1) return;
       
      if(sequence_ != space_.sequence()) 
      {
        buildMaps();

        // store actual sequence number 
        sequence_ = space_.sequence();
      }
      
      const int links = nlinks();
      // write buffers 
      for(int l=0; l<links; ++l) 
      {
        // reset buffers, keeps memory  
        buffer_[l].clear();

        writeBuffer( l , buffer_[l] , df.leakPointer());
      }

      // exchange data to other procs 
      buffer_ = mpAccess_.exchange( buffer_ );
      
      // read buffers 
      for(int l=0; l<links; ++l) 
      {
        readBuffer( l , buffer_[l] , df.leakPointer() , (OperationImp *) 0 );
      }
    }

  private:  
    // write data to object stream 
    template <class DataImp> 
    void writeBuffer(const int link,
                     ObjectStreamType & os, 
                     const DataImp* data) const 
    {
      const IndexMapType & indexMap = sendIndexMap_[ linkRank_ [link ] ]; 
      const int size = indexMap.size();

      // reserve buffer memory at once 
      os.reserve( size * sizeof(DataImp) );

      double val= 0.0;
      for(int i=0; i<size; ++i)
      {
        val = data[ indexMap[i] ];
        os.write( val );
      }
    }
  
    // read data from object stream to data vector 
    template <class DataImp, class OperationImp> 
    void readBuffer(const int link ,
                    ObjectStreamType & os, 
                    DataImp* data,
                    const OperationImp *) const 
    {
      const IndexMapType & indexMap = recvIndexMap_[ linkRank_ [link ] ]; 
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
    // type of communication manager object which does communication 
    typedef DependencyCache<SpaceImp> DependencyCacheType;

    typedef CommManagerSingletonKey<SpaceImp> KeyType;
    typedef CommManagerFactory<KeyType, DependencyCacheType> FactoryType;

    typedef SingletonList< KeyType , DependencyCacheType , FactoryType > CommunicationProviderType;

    typedef SpaceImp SpaceType;
    const SpaceType & space_; 

    const KeyType key_;

    const int mySize_;

    
    // is singleton per space 
    DependencyCacheType & cache_;
    CommunicationManager(const CommunicationManager & org);
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

    //! exchange discrete function to all procs we share data with 
    //! by using given OperationImp when receiving data from other procs 
    template <class DiscreteFunctionType> 
    void exchange(DiscreteFunctionType & df) 
    {
      cache_.exchange( df, (OperationImp*) 0 );
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
   
  //! \brief Default CommunicationManager class just using the grids communicate
  //! method 
  template <class SpaceImp, 
            class OperationImp = DFCommunicationOperation::Copy >
  class CommunicationManager 
  {
    typedef SpaceImp SpaceType; 
    typedef typename SpaceType :: GridPartType GridPartType; 

    // gridPart for communication 
    const GridPartType & gridPart_; 

    CommunicationManager(const CommunicationManager &);
  public:
    //! constructor taking space, but here only storing gridPart for
    //! communication
    CommunicationManager(const SpaceType & space)
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
  // end toggle AULGrid yes/no
#endif 
} // end namespace Dune 
#endif
