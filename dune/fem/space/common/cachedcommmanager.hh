#ifndef DUNE_FEM_CACHED_COMMUNICATION_MANAGER_HH
#define DUNE_FEM_CACHED_COMMUNICATION_MANAGER_HH

#include <cassert>
#include <cstddef>

//- system includes
#include <iostream>
#include <map>
#include <queue>
#include <memory>
#include <vector>

//- dune-common includes
#include <dune/common/math.hh>
#include <dune/common/timer.hh>
#include <dune/common/visibility.hh>

//- dune-grid includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/entitycommhelper.hh>

// include alugrid headers to have to communicator class from ALUGrid
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/3d/alu3dinclude.hh>
#endif

//- dune-fem includes
#include <dune/fem/common/hybrid.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/commindexmap.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class DiscreteFunctionSpace >
    class PetscDiscreteFunction;

    class IsDiscreteFunction;



    /** @addtogroup Communication Communication
        @{
     **/

// only if ALUGrid found and was build for parallel runs
// if HAVE_DUNE_ALUGRID is not defined, ALU3DGRID_PARALLEL shouldn't be either
#if ALU3DGRID_PARALLEL

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

    protected:
      // type of communication indices
      typedef CommunicationIndexMap IndexMapType;

      // type of IndexMapVector
      typedef IndexMapType *IndexMapVectorType;

      // type of set of links
      typedef std :: set< int > LinkStorageType;

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

      // ALUGrid communicatior Class
      MPAccessInterfaceType *mpAccess_;

      // exchange time
      double exchangeTime_;
      // setup time
      double buildTime_;

      //! know grid sequence number
      int sequence_;

      int nonBlockingObjects_ ;

    protected:
      template< class LinkStorage, class IndexMapVector, InterfaceType CommInterface >
      class LinkBuilder;

      /////////////////////////////////////////////////////////////////
      //  begin NonBlockingCommunication
      /////////////////////////////////////////////////////////////////

      class NonBlockingCommunication
      {
        typedef DependencyCache < Space > DependencyCacheType;

#if HAVE_DUNE_ALUGRID
        typedef MPAccessInterfaceType :: NonBlockingExchange NonBlockingExchange;

        template <class DiscreteFunction>
        class Pack : public NonBlockingExchange :: DataHandleIF
        {
        protected:
          NonBlockingCommunication& commObj_;
          const DiscreteFunction& discreteFunction_;

        public:
          Pack( NonBlockingCommunication& commObj, const DiscreteFunction& df )
          : commObj_( commObj ), discreteFunction_( df )
          {}

          void pack( const int link, ObjectStreamType& buffer )
          {
            commObj_.pack( link, buffer, discreteFunction_ );
          }

          void unpack( const int link, ObjectStreamType& buffer )
          {
            DUNE_THROW(InvalidStateException,"Pack::unpack should not be called!");
          }
        };

        template <class DiscreteFunction, class Operation>
        class Unpack : public NonBlockingExchange :: DataHandleIF
        {
        protected:
          NonBlockingCommunication& commObj_;
          DiscreteFunction& discreteFunction_;

          // communication operation (usually ADD or COPY)
          const Operation operation_;

        public:
          Unpack( NonBlockingCommunication& commObj, DiscreteFunction& df )
          : commObj_( commObj ), discreteFunction_( df ), operation_()
          {}

          void pack( const int link, ObjectStreamType& buffer )
          {
            DUNE_THROW(InvalidStateException,"Unpack::pack should not be called!");
          }

          void unpack( const int link, ObjectStreamType& buffer )
          {
            commObj_.unpack( link, buffer, discreteFunction_, operation_ );
          }
        };
#else   // ALUGRID_HAS_NONBLOCKING_COMM is false
        typedef int NonBlockingExchange;
#endif

        // create an unique tag for the communication
        DUNE_EXPORT static int getMessageTag()
        {
          enum { initial = 665 };
          static int tagCounter = initial ;
          ++ tagCounter;
          int messageTag = tagCounter ;

          // avoid overflow
          if( messageTag < 0 )
          {
            messageTag = initial ;
            tagCounter = initial ;
          }
          return messageTag;
        }

      public:
        NonBlockingCommunication( DependencyCacheType& dependencyCache, const int mySize )
          : dependencyCache_( dependencyCache ),
            nonBlockingExchange_(),
            buffer_(),
            exchangeTime_( 0.0 ),
            mySize_( mySize )
        {
          // update cache ( if necessary )
          dependencyCache_.rebuild();

          // notify dependency cache of open communication
          dependencyCache_.attachComm();
        }

        // copy constructor
        NonBlockingCommunication( const NonBlockingCommunication& other )
          : dependencyCache_( other.dependencyCache_ ),
            nonBlockingExchange_(),
            buffer_(),
            exchangeTime_( 0.0 ),
            mySize_( other.mySize_ )
        {
          // update cache ( if necessary )
          dependencyCache_.rebuild();

          // notify dependency cache of open communication
          dependencyCache_.attachComm();
        }

        ~NonBlockingCommunication()
        {
          // if this assertion fails some communication has not been finished
          assert( ! nonBlockingExchange_ );
          // notify dependency cache that comm is finished
          dependencyCache_.detachComm() ;
        }

        template < class DiscreteFunctionSpace >
        void send( const PetscDiscreteFunction< DiscreteFunctionSpace >& discreteFunction )
        {
          // nothing to do for the PetscDiscreteFunction here
        }

        template < class DiscreteFunction >
        void send( const DiscreteFunction& discreteFunction )
        {
          // check that object is in non-sent state
          assert( ! nonBlockingExchange_ );

          // on serial runs: do nothing
          if( mySize_ <= 1 ) return;

          // take time
          Dune::Timer sendTimer ;

          // this variable can change during rebuild
          const int nLinks = dependencyCache_.nlinks();

          // resize buffer vector
          buffer_.resize( nLinks );

#if HAVE_DUNE_ALUGRID
          // get non-blocking exchange object from mpAccess including message tag
          nonBlockingExchange_.reset( dependencyCache_.mpAccess().nonBlockingExchange( getMessageTag() ) );

          // pack data object
          Pack< DiscreteFunction > packData( *this, discreteFunction );

          // perform send operation including packing of data
          nonBlockingExchange_->send( buffer_, packData );
#else
          // write buffers
          for( int link = 0; link < nLinks; ++link )
            pack( link, buffer_[ link ], discreteFunction );
#endif

          // store time needed for sending
          exchangeTime_ = sendTimer.elapsed();
        }

        //! receive data for discrete function and given operation
        template < class DiscreteFunctionSpace, class Operation >
        double receive( PetscDiscreteFunction< DiscreteFunctionSpace >& discreteFunction,
                        const Operation& operation )
        {
          // take time
          Dune::Timer exchTimer;

          // PetscDiscreteFunction has it's own communication
          discreteFunction.dofVector().communicateNow( operation );

          return exchTimer.elapsed();
        }

        //! receive data for discrete function and given operation
        template < class DiscreteFunction, class Operation >
        double receive( DiscreteFunction& discreteFunction, const Operation& operation )
        {
          // on serial runs: do nothing
          if( mySize_ <= 1 ) return 0.0;

          // take time
          Dune::Timer recvTimer ;

#if HAVE_DUNE_ALUGRID
          // unpack data object
          Unpack< DiscreteFunction, Operation > unpackData( *this, discreteFunction );

          // receive data and unpack
          nonBlockingExchange_->receive( unpackData );
#else
          // use exchange for older ALUGrid versions (send and receive)
          buffer_ = dependencyCache_.mpAccess().exchange( buffer_ );

          // this variable can change during rebuild
          const int nLinks = buffer_.size();

          // read buffers and store to discrete function
          for( int link = 0; link < nLinks; ++link )
            unpack( link, buffer_[ link ], discreteFunction, operation );
#endif

          // store time needed for sending
          exchangeTime_ += recvTimer.elapsed();

#if HAVE_DUNE_ALUGRID
          // clear nonBlockingExchange object
          nonBlockingExchange_.reset();
#endif
          return exchangeTime_;
        }

        //! receive method with default operation
        template < class DiscreteFunction >
        double receive( DiscreteFunction& discreteFunction )
        {
          // get type of default operation
          typedef typename DiscreteFunction :: DiscreteFunctionSpaceType
            :: template CommDataHandle< DiscreteFunction > :: OperationType  DefaultOperationType;
          DefaultOperationType operation;
          return receive( discreteFunction, operation );
        }

      protected:
        template <class DiscreteFunction>
        void pack( const int link, ObjectStreamType& buffer, const DiscreteFunction& discreteFunction )
        {
          // reset buffer counters
          buffer.clear();
          // write data of discrete function to message buffer
          dependencyCache_.writeBuffer( link, buffer, discreteFunction );
        }

        template <class DiscreteFunction, class Operation>
        void unpack( const int link, ObjectStreamType& buffer,
                     DiscreteFunction& discreteFunction, const Operation& operation )
        {
          // read data of discrete function from message buffer
          dependencyCache_.readBuffer( link, buffer, discreteFunction, operation );
        }

      protected:
        DependencyCacheType& dependencyCache_;
        std::unique_ptr< NonBlockingExchange > nonBlockingExchange_ ;
        ObjectStreamVectorType buffer_;
        double exchangeTime_ ;
        const int mySize_;
      };

    public:
      typedef NonBlockingCommunication NonBlockingCommunicationType;

      //! return object for non-blocking communication
      NonBlockingCommunicationType nonBlockingCommunication()
      {
        // create non-blocking communication object
        return NonBlockingCommunicationType( *this, mySize_ );
      }
      /////////////////////////////////////////////////////////////////
      //  end NonBlockingCommunication
      /////////////////////////////////////////////////////////////////

      //! constructor taking space
      DependencyCache( const SpaceType &space, const InterfaceType interface, const CommunicationDirection dir )
      : space_( space ),
        gridPart_( space_.gridPart() ),
        interface_( interface ),
        dir_(dir),
        myRank_( gridPart_.comm().rank() ),
        mySize_( gridPart_.comm().size() ),
        linkStorage_(),
        recvIndexMap_( new IndexMapType[ mySize_ ] ),
        sendIndexMap_( new IndexMapType[ mySize_ ] ),
        mpAccess_( new MPAccessImplType( gridPart_.comm() ) ),
        exchangeTime_( 0.0 ),
        buildTime_( 0.0 ),
        sequence_( -1 ),
        nonBlockingObjects_( 0 )
      {}

      DependencyCache( const DependencyCache & ) = delete;

      //! destrcutor removeing mpAccess
      ~DependencyCache()
      {
        delete mpAccess_;
        mpAccess_ = 0;

        delete [] sendIndexMap_;
        sendIndexMap_ = 0;

        delete [] recvIndexMap_;
        recvIndexMap_ = 0;
      }

      //! return communication interface
      InterfaceType communicationInterface() const
      {
        return interface_;
      }

      //! return communcation direction
      CommunicationDirection communicationDirection() const
      {
        return dir_;
      }

      //! return time needed for last build
      double buildTime() const
      {
        return buildTime_;
      }

      //! return time needed for last exchange
      double exchangeTime() const
      {
        return exchangeTime_;
      }

      // build linkage and index maps
      inline void buildMaps();

      // notify for open non-blocking communications
      void attachComm()
      {
        ++ nonBlockingObjects_;
      }

      // notify for finished non-blocking communication
      void detachComm()
      {
        --nonBlockingObjects_;
        assert( nonBlockingObjects_ >= 0 );
      }

      bool noOpenCommunications() const
      {
        return true ;
      }
    protected:
      // check consistency of maps
      inline void checkConsistency();

      template< class LS, class IMV, InterfaceType CI >
      inline void buildMaps( LinkBuilder< LS, IMV, CI > &handle );

    public:
      //! return MPI rank of link
      inline int dest( const int link ) const
      {
        return mpAccess().dest()[ link ];
      }

      //! return number of links
      inline int nlinks() const
      {
        return mpAccess().nlinks();
      }

      //! check if grid has changed and rebuild cache if necessary
      inline void rebuild()
      {
        // only in parallel we have to do something
        if( mySize_ <= 1 ) return;

        // make sure all non-blocking communications have been finished by now
        assert( noOpenCommunications() );
#ifndef NDEBUG
        // make sure buildMaps is called on every process
        // otherwise the programs wait here until forever
        int willRebuild = (sequence_ != space_.sequence()) ? 1 : 0;
        const int myRebuild = willRebuild;

        // send willRebuild from rank 0 to all
        gridPart_.comm().broadcast( &willRebuild, 1 , 0);

        assert( willRebuild == myRebuild );
#endif

        // check whether grid has changed.
        if( sequence_ != space_.sequence() )
        {
          // take timer needed for rebuild
          Dune::Timer buildTime;

          buildMaps();
          sequence_ = space_.sequence();

          // store time needed
          buildTime_ = buildTime.elapsed();
        }
      }

      //! exchange data of discrete function
      template< class DiscreteFunction, class Operation >
      inline void exchange( DiscreteFunction &discreteFunction, const Operation& operation );

      //! write data of discrete function to buffer
      template< class DiscreteFunction >
      inline void writeBuffer( ObjectStreamVectorType &osv, const DiscreteFunction &discreteFunction ) const;

      //! read data of discrete function from buffer
      template< class DiscreteFunctionType, class Operation >
      inline void readBuffer( ObjectStreamVectorType &osv,
                              DiscreteFunctionType &discreteFunction,
                              const Operation& operation ) const;

      //! return reference to mpAccess object
      inline MPAccessInterfaceType &mpAccess()
      {
        assert( mpAccess_ );
        return *mpAccess_;
      }

      //! return reference to mpAccess object
      inline const MPAccessInterfaceType &mpAccess() const
      {
        assert( mpAccess_ );
        return *mpAccess_;
      }

    protected:
      // specialization for PetscDiscreteFunction doing nothing
      template< class DiscreteFunctionSpace >
      inline void writeBuffer( const int link,
                               ObjectStreamType &str,
                               const PetscDiscreteFunction< DiscreteFunctionSpace > &discreteFunction ) const
      {
        DUNE_THROW(NotImplemented,"writeBuffer not implemented for PetscDiscteteFunction" );
      }

      // write data of DataImp& vector to object stream
      // --writeBuffer
      template< class DiscreteFunction >
      inline void writeBuffer( const int link,
                               ObjectStreamType &str,
                               const DiscreteFunction &discreteFunction ) const
      {
        assert( sequence_ == space_.sequence() );
        const auto &indexMap = sendIndexMap_[ dest( link ) ];
        const int size = indexMap.size();

        typedef typename DiscreteFunction :: DofType DofType;

        // reserve write buffer for storage of dofs
        typename DiscreteFunction::DiscreteFunctionSpaceType::LocalBlockIndices localBlockIndices;
        str.reserve( size * Hybrid::size( localBlockIndices ) * sizeof( DofType ) );
        for( int i = 0; i < size; ++i )
        {
          const auto &block = discreteFunction.dofVector()[ indexMap[ i ] ];
          Hybrid::forEach( localBlockIndices, [ &str, &block ] ( auto &&k ) { str.writeUnchecked( block[ k ] ); } );
        }
      }

      // read data from object stream to DataImp& data vector
      // specialization for PetscDiscreteFunction doing nothing
      template< class DiscreteFunctionSpace, class Operation >
      inline void readBuffer( const int link,
                              ObjectStreamType &str,
                              PetscDiscreteFunction< DiscreteFunctionSpace > &discreteFunction,
                              const Operation& ) const
      {
        DUNE_THROW(NotImplemented,"readBuffer not implemented for PetscDiscteteFunction" );
      }

      // read data from object stream to DataImp& data vector
      // --readBuffer
      template< class DiscreteFunction, class Operation >
      inline void readBuffer( const int link,
                              ObjectStreamType &str,
                              DiscreteFunction &discreteFunction,
                              const Operation& operation ) const
      {
        static_assert( ! std::is_pointer< Operation > :: value,
                       "DependencyCache::readBuffer: Operation needs to be a reference!");

        assert( sequence_ == space_.sequence() );
        typedef typename DiscreteFunction :: DofType DofType;

        // get index map of rank belonging to link
        const auto &indexMap = recvIndexMap_[ dest( link ) ];

        const int size = indexMap.size();
        // make sure that the receive buffer has the correct size
        typename DiscreteFunction::DiscreteFunctionSpaceType::LocalBlockIndices localBlockIndices;
        assert( static_cast< std::size_t >( size * Hybrid::size( localBlockIndices ) * sizeof( DofType ) ) <= static_cast< std::size_t >( str.size() ) );
        for( int i = 0; i < size; ++i )
        {
          auto &&block = discreteFunction.dofVector()[ indexMap[ i ] ];
          Hybrid::forEach( localBlockIndices, [ &str, &operation, &block ] ( auto &&k ) {
              DofType value;
#if HAVE_DUNE_ALUGRID
              str.readUnchecked( value );
#else // #if HAVE_DUNE_ALUGRID
              str.read( value );
#endif // #else // #if HAVE_DUNE_ALUGRID
              // apply operation, i.e. COPY, ADD, etc.
              operation( value, block[ k ] );
            } );
        }
      }
    };

    // --LinkBuilder
    template< class Space >
    template< class LinkStorage, class IndexMapVector, InterfaceType CommInterface >
    class DependencyCache< Space > :: LinkBuilder
    : public CommDataHandleIF
      < LinkBuilder< LinkStorage, IndexMapVector, CommInterface >,
                     typename SpaceType :: BlockMapperType :: GlobalKeyType >
    {
    public:
      typedef LinkStorage LinkStorageType;

      typedef IndexMapVector IndexMapVectorType;

      typedef typename SpaceType :: BlockMapperType BlockMapperType;
      typedef typename BlockMapperType :: GlobalKeyType  GlobalKeyType;

      typedef GlobalKeyType DataType;

    protected:
      const GlobalKeyType myRank_;
      const GlobalKeyType mySize_;

      LinkStorageType &linkStorage_;

      IndexMapVectorType &sendIndexMap_;
      IndexMapVectorType &recvIndexMap_;

      const SpaceType &space_;
      const BlockMapperType &blockMapper_;

    public:
      LinkBuilder( LinkStorageType &linkStorage,
                   IndexMapVectorType &sendIdxMap,
                   IndexMapVectorType &recvIdxMap,
                   const SpaceType &space )
      : myRank_( space.gridPart().comm().rank() ),
        mySize_( space.gridPart().comm().size() ),
        linkStorage_( linkStorage ),
        sendIndexMap_( sendIdxMap ),
        recvIndexMap_( recvIdxMap ),
        space_( space ),
        blockMapper_( space.blockMapper() )
      {}

    protected:
      void sendBackSendMaps()
      {
        // create ALU communicator
        MPAccessImplType mpAccess( space_.gridPart().comm() );

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
          sendIndexMap_[ dest[link] ].writeToBuffer( osv[link] );

        // exchange data
        osv = mpAccess.exchange( osv );

        // read all send maps from buffer
        for(int link=0; link<nlinks; ++link)
          sendIndexMap_[ dest[link] ].readFromBuffer( osv[link] );
      }

    public:
      //! desctructor
      ~LinkBuilder()
      {
        sendBackSendMaps();
      }

      //! returns true if combination is contained
      bool contains( int dim, int codim ) const
      {
        return space_.blockMapper().contains( codim );
      }

      //! return whether we have a fixed size
      bool fixedSize( int dim, int codim ) const
      {
        return false;
      }

      //! read buffer and apply operation
      template< class MessageBuffer, class Entity >
      void gather( MessageBuffer &buffer, const Entity &entity ) const
      {
        // check whether we are a sending entity
        const auto myPartitionType = entity.partitionType();
        const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );

        // if we send data then send rank and dofs
        if( send )
        {
          // send rank for linkage
          buffer.write( myRank_ );

          const int numDofs = blockMapper_.numEntityDofs( entity );

          // int should be GlobalKey !!!!
          typedef std::vector< GlobalKeyType >  IndicesType ;
          IndicesType indices( numDofs );

          // copy all global keys
          blockMapper_.mapEachEntityDof( entity, AssignFunctor< IndicesType >( indices ) );

          // write global keys to message buffer
          for( int i = 0; i < numDofs; ++i )
            buffer.write( indices[ i ] );
        }
      }

      //! read buffer and apply operation
      template< class MessageBuffer, class Entity >
      void scatter( MessageBuffer &buffer, const Entity &entity, const size_t dataSize )
      {
        // if data size > 0 then other side is sender
        if( dataSize > 0 )
        {
          // read rank of other side
          GlobalKeyType rank;
          buffer.read( rank );
          assert( (rank >= 0) && (rank < mySize_) );

          // check whether we are a sending entity
          const auto myPartitionType = entity.partitionType();
          const bool receive = EntityCommHelper< CommInterface > :: receive( myPartitionType );

          // insert rank of link into set of links
          linkStorage_.insert( rank );

          // read indices from stream
          typedef std::vector< GlobalKeyType >  IndicesType ;
          IndicesType indices( dataSize - 1 );
          for(size_t i=0; i<dataSize-1; ++i)
            buffer.read( indices[i] );

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

            // map each entity dof and store in indices
            blockMapper_.mapEachEntityDof( entity, AssignFunctor< IndicesType >( indices ) );

            // insert receiving dofs
            recvIndexMap_[ rank ].insert( indices );
          }
        }
      }

      //! return local dof size to be communicated
      template< class Entity >
      size_t size( const Entity &entity ) const
      {
        const PartitionType myPartitionType = entity.partitionType();
        const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );
        return (send) ? (blockMapper_.numEntityDofs( entity ) + 1) : 0;
      }
    };



    template< class Space >
    inline void DependencyCache< Space > :: buildMaps()
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
        DUNE_THROW( NotImplemented, "DependencyCache for the given interface has not been implemented, yet." );
#ifndef NDEBUG
      // checks that sizes of index maps are equal on sending and receiving proc
      checkConsistency();
#endif
    }


    template< class Space >
    template< class LS, class IMV, InterfaceType CI >
    inline void DependencyCache< Space > :: buildMaps( LinkBuilder< LS, IMV, CI > &handle )
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
      mpAccess().insertRequestSymetric( linkStorage_ );
    }

    template< class Space >
    inline void DependencyCache< Space > :: checkConsistency()
    {
      const int nLinks = nlinks();

      ObjectStreamVectorType buffer( nLinks );

      // check that order and size are consistent
      for(int l=0; l<nLinks; ++l)
      {
        buffer[l].clear();
        const int sendSize = sendIndexMap_[ dest( l ) ].size();
        buffer[l].write( sendSize );
        for(int i=0; i<sendSize; ++i)
          buffer[l].write( i );
      }

      // exchange data to other procs
      buffer = mpAccess().exchange( buffer );

      // check that order and size are consistent
      for(int l=0; l<nLinks; ++l)
      {
        const int recvSize = recvIndexMap_[ dest( l ) ].size();
        int sendedSize;
        buffer[l].read( sendedSize );

        // compare sizes, must be the same
        if( recvSize != sendedSize )
        {
          DUNE_THROW(InvalidStateException,"Sizes do not match!" << sendedSize << " o|r " << recvSize);
        }

        for(int i=0; i<recvSize; ++i)
        {
          int idx;
          buffer[l].read( idx );

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
    inline void DependencyCache< Space > :: exchange( DiscreteFunction &discreteFunction, const Operation& operation )
    {
      // on serial runs: do nothing
      if( mySize_ <= 1 ) return;

      // create non-blocking communication object
      NonBlockingCommunicationType nbc( *this, mySize_ );

      // perform send operation
      nbc.send( discreteFunction );

      // store time for send and receive of data
      exchangeTime_ = nbc.receive( discreteFunction, operation );
    }

    template< class Space >
    template< class DiscreteFunction >
    inline void DependencyCache< Space > :: writeBuffer( ObjectStreamVectorType &osv,
                                                         const DiscreteFunction &discreteFunction ) const
    {
      const int numLinks = nlinks();
      for( int link = 0; link < numLinks; ++link )
        writeBuffer( link, osv[ link ], discreteFunction );
    }

    template< class Space >
    template< class DiscreteFunction, class Operation >
    inline void DependencyCache< Space > :: readBuffer( ObjectStreamVectorType &osv, DiscreteFunction &discreteFunction,
                                                        const Operation& operation ) const
    {
      const int numLinks = nlinks();
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
        : space_(space), interface_(interface), dir_(dir)
      {}

      //! copy constructor
      CommManagerSingletonKey(const CommManagerSingletonKey & org)
        : space_(org.space_), interface_(org.interface_), dir_(org.dir_)
      {}

      //! returns true if indexSet pointer and numDofs are equal
      bool operator == (const CommManagerSingletonKey & otherKey) const
      {
        // mapper of space is singleton
        return (&(space_.blockMapper()) == & (otherKey.space_.blockMapper()) );
      }

      //! return reference to index set
      const SpaceImp & space() const
      {
        return space_;
      }

      //! return communication interface
      InterfaceType interface() const
      {
        return interface_;
      }

      //! return communication direction
      CommunicationDirection direction() const
      {
        return dir_;
      }
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
        return new ObjectImp(key.space(), key.interface(), key.direction());
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
      DependencyCacheType &cache_;
      CommunicationManager(const ThisType& org);
    public:
      // type of non-blocking communication object
      typedef typename DependencyCacheType :: NonBlockingCommunicationType  NonBlockingCommunicationType;

      //! constructor taking space and communication interface/direction
      CommunicationManager(const SpaceType & space,
                           const InterfaceType interface,
                           const CommunicationDirection dir)
        : space_(space)
        , key_(space_,interface,dir)
        , mySize_(space_.gridPart().comm().size())
        , cache_(CommunicationProviderType::getObject(key_))
      {}

      //! constructor taking space and communication interface/direction
      CommunicationManager(const SpaceType & space)
        : space_(space)
        , key_(space_, space.communicationInterface(), space.communicationDirection())
        , mySize_(space_.gridPart().comm().size())
        , cache_(CommunicationProviderType::getObject(key_))
      {}

      //! remove object comm
      ~CommunicationManager()
      {
        CommunicationProviderType::removeObject(cache_);
      }

      //! return communication interface
      InterfaceType communicationInterface() const
      {
        return cache_.communicationInterface();
      }

      //! return communcation direction
      CommunicationDirection communicationDirection() const
      {
        return cache_.communicationDirection();
      }

      //! return time needed for last build
      double buildTime() const
      {
        return cache_.buildTime();
      }

      //! return time needed for last exchange
      double exchangeTime() const
      {
        return cache_.exchangeTime();
      }

      MPAccessInterfaceType& mpAccess()
      {
        return cache_.mpAccess();
      }

      //! return object for non-blocking communication
      NonBlockingCommunicationType nonBlockingCommunication() const
      {
        return cache_.nonBlockingCommunication();
      }

      //! exchange discrete function to all procs we share data
      //! using the copy operation
      template <class DiscreteFunctionType>
      void exchange(DiscreteFunctionType & df) const
      {
        // get type of default operation
        typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
          :: template CommDataHandle< DiscreteFunctionType > :: OperationType  DefaultOperationType;

        // create default operation
        DefaultOperationType operation;

        exchange( df, operation );
      }

      //! exchange discrete function to all procs we share data
      //! using the given operation
      template <class DiscreteFunctionType, class Operation>
      void exchange(DiscreteFunctionType & df, const Operation& operation ) const
      {
        cache_.exchange( df, operation );
      }

      //! write given df to given buffer
      template <class ObjectStreamVectorType, class DiscreteFunctionType>
      void writeBuffer(ObjectStreamVectorType& osv, const DiscreteFunctionType & df) const
      {
        cache_.writeBuffer( osv, df );
      }

      // read given df from given buffer
      template <class ObjectStreamVectorType, class DiscreteFunctionType>
      void readBuffer(ObjectStreamVectorType& osv, DiscreteFunctionType & df) const
      {
        typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
          :: template CommDataHandle<DiscreteFunctionType> :: OperationType OperationType;

        // communication operation to be performed on the received data
        OperationType operation;

        readBuffer( osv, df, operation );
      }

      // read given df from given buffer
      template <class ObjectStreamVectorType, class DiscreteFunctionType, class OperationType>
      void readBuffer(ObjectStreamVectorType& osv, DiscreteFunctionType & df, const OperationType& operation) const
      {
        cache_.readBuffer( osv, df , operation);
      }

      //! rebuild underlying cache if necessary
      void rebuildCache()
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
        DiscreteFunctionCommunicatorInterface()
        {}
      public:
        virtual ~DiscreteFunctionCommunicatorInterface()
        {}

        virtual MPAccessType& mpAccess() = 0;
        virtual void writeBuffer(ObjectStreamVectorType&) const = 0;
        virtual void readBuffer(ObjectStreamVectorType&) = 0;
        virtual void rebuildCache() = 0;

        virtual bool handles ( IsDiscreteFunction &df ) const = 0;
      };

      //! communicated object implementation
      //! default operation is copy, because these lists are used to
      //! restore consistency only
      template <class DiscreteFunctionImp,
                class MPAccessType,
                class ObjectStreamVectorType,
                class OperationType >
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

        //! operation performed on received data
        const OperationType operation_;
      public:
        //! constructor taking disctete function
        DiscreteFunctionCommunicator(DiscreteFunctionType& df, const OperationType& op )
          : df_(df), comm_(df_.space()), operation_( op )
        {}

        //! return ALUGrid communicator
        virtual MPAccessType& mpAccess()
        {
          return comm_.mpAccess();
        }

        //! write discrete function to all buffers
        virtual void writeBuffer(ObjectStreamVectorType& osv) const
        {
          comm_.writeBuffer(osv,df_);
        }

        //! read discrete function from all buffers
        virtual void readBuffer(ObjectStreamVectorType& osv)
        {
          comm_.readBuffer(osv, df_, operation_ );
        }

        //! rebuild cache if grid changed
        virtual void rebuildCache()
        {
          comm_.rebuildCache();
        }

        virtual bool handles ( IsDiscreteFunction &df ) const { return (&static_cast< IsDiscreteFunction & >( df_ ) == &df); }
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
      typedef std::list < std::unique_ptr< CommObjInterfaceType > > CommObjListType;
      CommObjListType objList_;

      // number of processors
      int mySize_;

    public:
      //! constructor creating list of communicated objects
      template <class CombinedObjectType>
      CommunicationManagerList(CombinedObjectType& cObj) :
        mySize_( -1 )
      {
        // add all discrete functions containd in cObj to list
        cObj.addToList(*this);
      }

      //! default constructor
      CommunicationManagerList()
        : mySize_( -1 )
      {}

      CommunicationManagerList ( const CommunicationManagerList & ) = delete;

      //! add one discrete function to the list with given unpack operation
      template <class DiscreteFunctionImp, class Operation>
      void addToList(DiscreteFunctionImp &df, const Operation& operation )
      {
        // type of communication object
        typedef DiscreteFunctionCommunicator<DiscreteFunctionImp,
                                             MPAccessInterfaceType,
                                             ObjectStreamVectorType,
                                             Operation > CommObj;
        CommObj * obj = new CommObj(df, operation);
        objList_.push_back( std::unique_ptr< CommObjInterfaceType > (obj) );

        // if mySize wasn't set, set to number of processors
        if( mySize_ < 0 )
        {
          // get ALUGrid communicator
          MPAccessInterfaceType& mpAccess = objList_.front()->mpAccess();

          // set number of processors
          mySize_ = mpAccess.psize();
        }
      }

      //! add one discrete function to the list
      template <class DiscreteFunctionImp>
      void addToList(DiscreteFunctionImp &df)
      {
        DFCommunicationOperation::Copy operation;
        addToList( df, operation );
      }

      template< class DiscreteFunction >
      void removeFromList ( DiscreteFunction &df )
      {
        const auto handles = [ &df ] ( const std::unique_ptr< CommObjInterfaceType > &commObj ) { assert( commObj ); return commObj->handles( df ); };
        CommObjListType::reverse_iterator pos = std::find_if( objList_.rbegin(), objList_.rend(), handles );
        if( pos != objList_.rend() )
          objList_.erase( --pos.base() );
        else
          DUNE_THROW( RangeError, "Trying to remove discrete function that was never added" );
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
          // rebuild cahce if grid has changed
          for(auto& elem : objList_)
            elem->rebuildCache();

          // get ALUGrid communicator
          auto& mpAccess = objList_.front()->mpAccess();

          // create buffer
          ObjectStreamVectorType osv( mpAccess.nlinks() );

          // fill buffers
          for(auto& elem : objList_)
            elem->writeBuffer(osv);

          // exchange data
          osv = mpAccess.exchange(osv);

          // read buffers
          for(auto& elem : objList_)
            elem->readBuffer(osv);
        }
      }
    };
#endif  // #if ALU3DGRID_PARALLEL
    //@}

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_CACHED_COMMUNICATION_MANAGER_HH
