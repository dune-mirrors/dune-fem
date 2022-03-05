#ifndef DUNE_FEM_COMMUNICATION_MANAGER_HH
#define DUNE_FEM_COMMUNICATION_MANAGER_HH

#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/storage/singletonlist.hh>

// include ALUGrid to check whether the
// parallel version is avaiable
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/3d/alugrid.hh>
#endif

// default is: enabled
#ifndef WANT_CACHED_COMM_MANAGER
#define WANT_CACHED_COMM_MANAGER 1
#endif

#if ALU3DGRID_PARALLEL && WANT_CACHED_COMM_MANAGER
#define USE_CACHED_COMM_MANAGER
#else
#ifndef NDEBUG
#if HAVE_MPI == 0
  #ifdef DUNE_DEVEL_MODE
    #warning "HAVE_MPI == 0, therefore default CommunicationManager is used!"
  #endif
#elif !ALU3DGRID_PARALLEL
  #warning "No Parallel ALUGrid found, using default CommunicationManager!"
#elif ! WANT_CACHED_COMM_MANAGER
  #warning "CachedCommunication Manager disabled by WANT_CACHED_COMM_MANAGER=0!"
#endif
#endif
#endif

#undef WANT_CACHED_COMM_MANAGER

#ifdef USE_CACHED_COMM_MANAGER
#include "cachedcommmanager.hh"
#endif

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class DiscreteFunctionSpace >
    class PetscDiscreteFunction;

    class IsDiscreteFunction;

    class IsBlockVector;


    /** @addtogroup Communication Communication
        @{
    **/

    /** \class DefaultCommunicationManager
     *  \ingroup Communication
     *  \brief default communication manager using just the grids communicate
     *         method
     */
    template< class Space >
    class DefaultCommunicationManager
    {
    public:
      typedef Space SpaceType;

    protected:
      typedef DefaultCommunicationManager< Space > ThisType;

      /////////////////////////////////////////////////////////////////
      //  begin NonBlockingCommunication
      /////////////////////////////////////////////////////////////////
      class NonBlockingCommunication
      {
        const SpaceType& space_;
        const InterfaceType interface_;
        const CommunicationDirection dir_;

      public:
        NonBlockingCommunication( const SpaceType& space,
                                  InterfaceType interface,
                                  CommunicationDirection dir )
          : space_( space ),
            interface_( interface ),
            dir_ ( dir )
        {}

        //! send data for given discrete function
        template < class DiscreteFunction >
        void send( const DiscreteFunction& discreteFunction )
        {
          // nothing to do here, since DUNE does not support
          // non-blocking communication yet
        }

        //! receive data for discrete function and given operation
        template < class DiscreteFunctionSpace, class Operation >
        double receive( PetscDiscreteFunction< DiscreteFunctionSpace > & discreteFunction,
                        const Operation& operation )
        {
          // on serial runs: do nothing
          if( space_.gridPart().comm().size() <= 1 )
            return 0.0;

          // get stopwatch
          Dune::Timer exchangeT;

          // PetscDiscreteFunction has it's own communication
          discreteFunction.dofVector().communicateNow( operation );

          return exchangeT.elapsed();
        }

        //! receive data for discrete function and given operation
        template < class DiscreteFunction, class Operation >
        double receive( DiscreteFunction& discreteFunction, const Operation& operation )
        {
          // on serial runs: do nothing
          if( space_.gridPart().comm().size() <= 1 )
            return 0.0;

          // get type of data handle from the discrete function space
          typedef typename DiscreteFunction
            :: template CommDataHandle< Operation > :: Type
            DataHandleType;

          // get stopwatch
          Dune::Timer exchangeT;

          // communicate data
          DataHandleType dataHandle = discreteFunction.dataHandle( operation );
          space_.gridPart().communicate( dataHandle, interface_ , dir_ );

          // store time
          return exchangeT.elapsed();
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

      };

      const SpaceType& space_;

      const InterfaceType interface_;
      const CommunicationDirection dir_;

      mutable double exchangeTime_;

    public:
      typedef NonBlockingCommunication  NonBlockingCommunicationType;

      //! constructor taking space and communication interface/direction
      explicit DefaultCommunicationManager
        ( const SpaceType &space,
          const InterfaceType interface,
          const CommunicationDirection dir)
      : space_( space ),
        interface_( interface ),
        dir_ ( dir ),
        exchangeTime_(0.0)
      {}

      DefaultCommunicationManager ( const DefaultCommunicationManager& ) = delete;

      /** \brief return communication interface */
      InterfaceType communicationInterface() const {
        return interface_;
      }

      /** \brief return communication direction */
      CommunicationDirection communicationDirection() const
      {
        return dir_;
      }

      /** \brief return time needed for last build

          \return time needed for last build of caches (if needed)
      */
      double buildTime() const { return 0.0; }

      /** \brief return time needed for last exchange of data

          \return time needed for last exchange of data
      */
      double exchangeTime() const { return exchangeTime_; }

      /** \brief return object for non-blocking communication

          \return NonBlockingCommunicationType containing send and receive facilities
        */
      NonBlockingCommunicationType nonBlockingCommunication() const
      {
        return NonBlockingCommunicationType( space_, interface_, dir_ );
      }

      /** \brief exchange data for a discrete function using the copy operation
       *
       *  \param  discreteFunction  discrete function to communicate
       */
      template< class DiscreteFunction >
      inline void exchange ( DiscreteFunction &discreteFunction ) const
      {
        // get type of default operation
        typedef typename DiscreteFunction :: DiscreteFunctionSpaceType ::
          template CommDataHandle< DiscreteFunction > :: OperationType DefaultOperationType;

        DefaultOperationType operation;
        exchange( discreteFunction, operation );
      }

      /** \brief exchange data for a discrete function using the given operation
       *
       *  The used operation is derived from the type of the op-pointer. The
       *  actual pointer is not used.
       *
       *  \param      discreteFunction  discrete function to communicate
       *  \param[in]  operation         a (phony) pointer to an operation
       */
      template< class DiscreteFunction, class Operation >
      inline void exchange ( DiscreteFunction &discreteFunction,
                             const Operation &operation ) const
      {
        // on serial runs: do nothing
        if( space_.gridPart().comm().size() <= 1 )
          return;

        NonBlockingCommunicationType nbc( space_, interface_, dir_ );

        // send data (probably done in receive)
        nbc.send( discreteFunction );

        exchangeTime_ = nbc.receive( discreteFunction, operation );
      }

      /** \brief exchange data for a discrete function using the given operation
       *
       *  The used operation is derived from the type of the op-pointer. The
       *  actual pointer is not used.
       *
       *  \param      discreteFunction  discrete function to communicate
       *  \param[in]  operation         a (phony) pointer to an operation
       */
      template< class DiscreteFunction, class Operation >
      inline void exchange ( const Space& space,
                             DiscreteFunction &discreteFunction,
                             const Operation &operation ) const
      {
#ifndef USE_CACHED_COMM_MANAGER
        DUNE_THROW(NotImplemented,"Communicating BlockVectorInterface and derived only works with cached communication!");
#endif
      }
    };



#ifndef USE_CACHED_COMM_MANAGER
    // if no ALUGrid found, supply default implementation
    //! \brief use Default CommunicationManager as Communication Manager
    template <class SpaceImp>
    class CommunicationManager
    : public DefaultCommunicationManager<SpaceImp>
    {
      typedef DefaultCommunicationManager<SpaceImp> BaseType;
      CommunicationManager(const CommunicationManager &);
    public:
      //! constructor taking space and communication interface/direction
      CommunicationManager(const SpaceImp & space,
          const InterfaceType interface,
          const CommunicationDirection dir)
        : BaseType(space,interface,dir)
      {}
      //! constructor taking space
      CommunicationManager(const SpaceImp & space)
        : BaseType(space,
                   space.communicationInterface(),
                   space.communicationDirection() )
      {}
    };



    //! Proxy class to DependencyCache which is singleton per space
    class CommunicationManagerList
    {
      //! communicated object interface
      class DiscreteFunctionCommunicatorInterface
      {
      protected:
        DiscreteFunctionCommunicatorInterface () = default;
      public:
        virtual ~DiscreteFunctionCommunicatorInterface () = default;
        virtual void exchange () const = 0;
        virtual bool handles ( IsDiscreteFunction &df ) const = 0;
      };

      //! communicated object implementation
      template <class DiscreteFunctionImp, class Operation>
      class DiscreteFunctionCommunicator
      : public DiscreteFunctionCommunicatorInterface
      {
        typedef DiscreteFunctionImp DiscreteFunctionType;
        typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType;

        DiscreteFunctionType& df_;
        CommunicationManagerType comm_;
        const Operation& operation_;

      public:
        //! constructor taking disctete function
        DiscreteFunctionCommunicator(DiscreteFunctionType& df, const Operation& op)
          : df_(df), comm_(df_.space()), operation_( op )
        {
        }

        // exchange discrete function
        void exchange () const
        {
          comm_.exchange( df_, operation_ );
        }

        bool handles ( IsDiscreteFunction &df ) const { return (&df_ == &df); }
      };

      typedef DiscreteFunctionCommunicatorInterface CommObjIFType;
      typedef std::list < std::unique_ptr< DiscreteFunctionCommunicatorInterface > > CommObjListType;
      CommObjListType objList_;

      CommunicationManagerList(const CommunicationManagerList&);
    public:
      CommunicationManagerList () = default;

      //! constructor
      template <class CombinedObjectType>
      CommunicationManagerList(CombinedObjectType& cObj)
      {
        cObj.addToList(*this);
      }

      //! add discrete function to communication list
      template <class DiscreteFunctionImp, class Operation>
      void addToList(DiscreteFunctionImp &df, const Operation& operation )
      {
        typedef DiscreteFunctionCommunicator<DiscreteFunctionImp, Operation> CommObjType;
        CommObjType* obj = new CommObjType( df, operation );
        objList_.push_back( std::unique_ptr< DiscreteFunctionCommunicatorInterface> (obj) );
      }

      //! add discrete function to communication list
      template <class DiscreteFunctionImp>
      void addToList(DiscreteFunctionImp &df)
      {
        DFCommunicationOperation::Copy operation;
        addToList( df, operation );
      }

      template< class DiscreteFunction >
      void removeFromList ( DiscreteFunction &df )
      {
        const auto handles = [ &df ] ( const std::unique_ptr< DiscreteFunctionCommunicatorInterface > &commObj ) { return commObj->handles( df ); };
        CommObjListType::reverse_iterator pos = std::find_if( objList_.rbegin(), objList_.rend(), handles );
        if( pos != objList_.rend() )
          objList_.erase( --pos.base() );
        else
          DUNE_THROW( RangeError, "Trying to remove discrete function that was never added" );
      }

      //! exchange discrete function to all procs we share data with
      //! by using given OperationImp when receiving data from other procs
      void exchange() const
      {
        typedef CommObjListType :: const_iterator iterator;
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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMUNICATION_MANAGER_HH
