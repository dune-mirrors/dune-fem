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
#include <dune/fem/discretefunction/common/dfcommunication.hh>

namespace Dune { 

  // only if ALUGrid found and was build for parallel runs 
#if HAVE_ALUGRID && ALU3DGRID_PARALLEL
  template <class SpaceImp, 
            class OperationImp = DFCommunicationOperation::Copy >
  class CommunicationManager 
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
        const bool interiorEn = en.partitionType() == InteriorEntity; 
        
        // read links and insert mapping 
        DataType val;
        
        buff.read( val );  
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

    typedef SpaceImp SpaceType; 
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

    // ALUGrid send/recv buffers 
    typedef ALU3DSPACE ObjectStream ObjectStreamType; 
    // ALUGrid communicatior Class 
    ALU3DSPACE MpAccessMPI mpAccess_;

    int sequence_;
    std::vector < int > linkRank_;
    std::vector< ObjectStreamType > osv_;
    // do not copy this class 
    CommunicationManager(const CommunicationManager &);
  public:
    //! constructor taking space 
    CommunicationManager(const SpaceType & space)
      : space_(space) , gridPart_(space_.gridPart()) 
      , myRank_(gridPart_.grid().comm().rank())
      , mySize_(gridPart_.grid().comm().size())
      , linkStorage_()
      , recvIndexMap_(mySize_)
      , sendIndexMap_(mySize_)
      , mpAccess_(MPI_COMM_WORLD)
      , sequence_(-1)
      , linkRank_()
    {
    }

  //! exchange discrete function to all procs we share data with 
  //! by using given OperationImp when receiving data from other procs 
  template <class DiscreteFunctionType> 
  void exchange(DiscreteFunctionType & df) 
  {
    if(sequence_ != space_.sequence()) buildMaps();
    
    const int links = mpAccess_.nlinks();
    // write buffers 
    for(int l=0; l<links; ++l) 
    {
      writeBuffer(osv_[l], sendIndexMap_[ linkRank_[l] ], df.leakPointer());
    }

    // exchange data to other procs 
    osv_ = mpAccess_.exchange( osv_ );
    
    // read buffers 
    for(int l=0; l<links; ++l) 
    {
      readBuffer(osv_[l], recvIndexMap_[ linkRank_[l] ], df.leakPointer() );
    }
  }

  private: 
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
        gridPart_.communicate( handle, All_All_Interface , BackwardCommunication );
      }

      // remove old linkage 
      mpAccess_.removeLinkage(); 
      // create new linkage 
      mpAccess_.insertRequestSymetric ( linkStorage_ );

      // get real rank numbers for each link 
      linkRank_ = mpAccess_.dest();

      // resize buffer vector 
      osv_.resize( mpAccess_.nlinks() );

      // store actual sequence number 
      sequence_ = space_.sequence();
    }

    // write data to object stream 
    void writeBuffer(ObjectStreamType & os, 
                     const IndexMapType & indexMap , 
                     const double * data) const 
    {
      const int size = indexMap.size();
      for(int i=0; i<size; ++i)
      {
        os.write( data[ indexMap[i] ] );
      }
    }
  
    // read data from object stream to data vector 
    void readBuffer(ObjectStreamType & os, 
                    const IndexMapType & indexMap , 
                    double * data) const 
    {
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
#else 
#warning "No Parallel ALUGrid found, using default CommunicationManager!"
  // if no ALUGrid found, supply default implementation 
   
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
