#ifndef DUNE_FEM_PASS_COMMON_PASS_HH
#define DUNE_FEM_PASS_COMMON_PASS_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/pass.hh> instead!"
#endif

#include <limits>
#include <string>
#include <tuple>
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/visibility.hh>

#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {

    // empty non-blocking communication for start pass as default argument
    struct EmptyNonBlockingComm
    {
      // initialize communication
      template <class Destination>
      void initComm( const Destination& ) const {}

      // receive data of previously initialized communication
      template <class Destination>
      void receiveComm( const Destination& ) const {}

      // cleanup overwritten data, i.e. ghost values
      template <class Destination>
      void finalizeComm( const Destination& ) const {}
    };

    /*! @addtogroup Pass
     *
     */
    /*!
     * @brief End marker for a compile-time list of passes.
     *
     */
    template < class ArgumentImp , int passIdImp,
               class NonBlockingCommunication = EmptyNonBlockingComm  >
    class StartPass
    {
      NonBlockingCommunication nonBlockingComm_;
    public:
      //! position in pass tree (0 for start pass)
      static const int passNum = 0;
      static const int passId = passIdImp;

      //! pass ids up to here (tuple of integral constants)
      typedef std::tuple< std::integral_constant< int, passIdImp > > PassIds;

      //! The argument (and destination) type of the overall operator
      typedef ArgumentImp GlobalArgumentType;
      //! End marker for tuple of return types of the passes
      typedef std::tuple<> NextArgumentType;

    public:
      //! empty constructor
      StartPass() : nonBlockingComm_() {}
      //! copy constructor
      StartPass(const StartPass&) : nonBlockingComm_() {}

      //- Public methods
      //! The pass method initialized the communication only
      void pass( const GlobalArgumentType& arg ) const
      {
        nonBlockingComm_.initComm( arg );
      }

      //! receive data for previously initialized communication
      void receiveCommunication( const GlobalArgumentType& arg ) const
      {
        nonBlockingComm_.receiveComm( arg );
      }

      //! cleanup of overwritten data. i.e. ghost values if neccessary
      void finalizeCommunication( const GlobalArgumentType& arg ) const
      {
        nonBlockingComm_.finalizeComm( arg );
      }

      //! Returns the closure of the destination tuple.
      NextArgumentType localArgument() const { return NextArgumentType(); }

      //! No memory needs to be allocated.
      void allocateLocalMemory() {}

      //! StartPass does not need a time
      void setTime(const double) {}

      //! return time step estimate
      double timeStepEstimate() const
      {
        return std::numeric_limits<double>::max();
      }
    };

    /**
     * @brief Base class for specific pass implementations.
     *
     * Pass not only provides the interface for the specialised implementations,
     * but also organizes the calls to other passes and assembles the results
     * of the preceding passes in a tuple of discrete functions which can be
     * used in the computations of the pass. The computations must be implemented
     * in the compute method of the derived classes.
     */
    template <class DiscreteModelImp, class PreviousPassImp , int passIdImp>
    class Pass :
      public Operator<typename PreviousPassImp::GlobalArgumentType,
                      typename DiscreteModelImp::Traits::DestinationType>
    {
      template <class PT, class PP, int PI>
      friend class Pass;
    public:
      //! little interface class for deleting discrete function
      //! held by this class
      template <class ObjectToDelete>
      class DeleteHandler
      {
      protected:
        // don't create intances of this class
        DeleteHandler () {}
      public:
        //! destructor
        virtual ~DeleteHandler () {}
        //! default implementation just deletes obj
        virtual void freeLocalMemory(ObjectToDelete * obj)
        {
          delete obj;
        }

        //! return reference to default object deleter
        DUNE_EXPORT static DeleteHandler<ObjectToDelete>& instance ()
        {
          static DeleteHandler<ObjectToDelete> mh;
          return mh;
        }
      };

      //! Type of the preceding pass.
      typedef PreviousPassImp PreviousPassType;

      //! position in pass tree
      static const int passNum = PreviousPassType::passNum + 1;
      static const int passId  = passIdImp;

      //! pass ids up to here (tuple of integral constants)
      typedef typename Dune::PushBackTuple< typename PreviousPassType::PassIds, std::integral_constant< int, passIdImp > >::type PassIds;

      //! Type of the discrete function which stores the result of this pass'
      //! computations.
      typedef typename DiscreteModelImp::Traits::DestinationType DestinationType;

      //! type of mem handler, which deletes destination
      typedef DeleteHandler<DestinationType> DeleteHandlerType;

      typedef typename DestinationType :: DiscreteFunctionSpaceType :: CommunicationManagerType
            :: NonBlockingCommunicationType  NonBlockingCommunicationType;

      //! Type of the discrete function which is passed to the overall operator
      //! by the user
      typedef typename PreviousPassType::GlobalArgumentType GlobalArgumentType;
      //! Tuple containing destination types of all preceding passes.
      typedef typename PreviousPassType::NextArgumentType LocalArgumentType;
      //! Tuple containing destination types of all preceding passes plus the
      //! global argument. This serves as the argument for this pass'
      //! computations
      typedef typename PushFrontTuple< LocalArgumentType, const GlobalArgumentType* >::type TotalArgumentType;
      //! Tuple containing destination types of all passes up to this one.
      typedef typename PushBackTuple< LocalArgumentType, DestinationType* >::type NextArgumentType;

    public:
      // return pass number
      int passNumber() const { return passNum; }

      //- Public methods
      //! Constructor
      //! \param pass Previous pass
      Pass(PreviousPassType& pass) :
        destination_(0),
        deleteHandler_(0),
        previousPass_(pass),
        time_(0.0),
        finalizeCommunication_( true )
      {
        // this ensures that the last pass doesn't allocate temporary memory
        // (its destination discrete function is provided from outside).
        previousPass_.allocateLocalMemory();
      }

      //! Destructor
      virtual ~Pass()
      {
        // if deleteHandler was set by derived class,
        // then use to delete destination_
        if( deleteHandler_ ) deleteHandler_->freeLocalMemory(destination_);
        destination_ = 0;
      }

      //! printTex info of operator
      void printTexInfo(std::ostream& out) const
      {
        previousPass_.printTexInfo(out);
      }

      //! \brief Application operator.
      //! The application operator is called by the client directly. It makes
      //! only sense to call this operator directly on the last pass.
      void operator()(const GlobalArgumentType& arg, DestinationType& dest) const
      {
        previousPass_.pass(arg);
        LocalArgumentType prevArg = previousPass_.localArgument();
        const TotalArgumentType totalArg = tuple_push_front( prevArg, &arg );
        this->compute(totalArg, dest);

        // if initComm has not been called for this pass, we have to
        // call finalizeCommunication for all previous passes
        if( finalizeCommunication_ )
          finalizeCommunication( arg );
      }
      //! Allocates the local memory of a pass, if needed.
      //! If memory is allocated, then deleteHandler must be set for removal of
      //! memory to avoid leaks
      virtual void allocateLocalMemory() = 0;

      //! Set time provider (which gives you access to the global time).
      void setTime(const double t)
      {
        previousPass_.setTime(t);
        time_ = t;
      }

      /** \brief return time step estimate for explicit Runge Kutta solver, calls
           recursively the method timeStepEstimateImpl of all previous passes.
           Make sure to overload the method timeStepEstimateImpl in your
           implementation if this method really does something. */
      double timeStepEstimate() const
      {
        double ret= std::min( previousPass_.timeStepEstimate(),
                              this->timeStepEstimateImpl() );
        return ret;
      }

      //! return current time of calculation
      double time() const { return time_; }

      //! return reference to internal discrete function
      const DestinationType & destination () const
      {
        assert(destination_);
        return *destination_;
      }

    public:
      //! Same as application operator, but uses own memory instead of the
      //! discrete function provided by the client. This method is called on all
      //! passes except the last one.
      void pass(const GlobalArgumentType& arg) const
      {
        // send my destination data needed by the next pass
        initComm();
        // since initComm was called we are not the last pass
        // and thus must not call finalizeCommunication
        finalizeCommunication_ = false ;
        operator()(arg, *destination_);
      }

      //! Returns a compilation of the results of the preceding passes
      NextArgumentType localArgument () const
      {
        typename PreviousPassType::NextArgumentType nextArg( previousPass_.localArgument() );
        return tuple_push_back( nextArg, destination_ );
      }

      /** \brief finalizeCommunication collects possbily initiated non-blocking
                 communications for all passes including the global argument
                 this method will be called from the next pass
      */
      void finalizeCommunication(const GlobalArgumentType& arg) const
      {
        // we only need the first argument
        // the other argument are known to the previous passes
        previousPass_.finalizeCommunication( arg );

        // this method has to be overloaded in the pass implementation
        finalizeComm();
        // reset finalizeCommunication flag
        finalizeCommunication_ = true;
      }

      /** \brief finalizeCommunication collects possbily initiated non-blocking
                 communications for all passes including the global argument
                 this method will be called from the next pass
      */
      void receiveCommunication(const GlobalArgumentType& arg) const
      {
        // we only need the first argument
        // the other argument are known to the previous passes
        previousPass_.receiveCommunication( arg );

        // this method has to be overloaded in the pass implementation
        receiveComm();
      }

      //! derived passes have to implement this method
      //! returning the time step estimate
      virtual double timeStepEstimateImpl() const
      {
        return std::numeric_limits<double>::max();
      }

      /** \brief requireCommunication returns true if the pass needs communication at all
       *         \note The default implementation returns \b true \b
       */
      virtual bool requireCommunication () const { return true; }

    protected:
      //! Does the actual computations. Needs to be overridden in the derived
      //! clases
      virtual void compute(const TotalArgumentType& arg,
                           DestinationType& dest) const = 0;

      /** \brief finalizeCommunication collects possbily initiated non-blocking
                 communications for all passes
      */
      void finalizeCommunication( const TotalArgumentType& totalArg ) const
      {
        // this method is called on the last pass which needs no
        // finalizing of communication, so call the correct method
        // on the previous pass
        // std::get<0> ( totalArg ) is the global argument type
        previousPass_.finalizeCommunication( *(std::get<0> ( totalArg )) );
      }

      /** \brief receiveCommunication collects possbily initiated non-blocking
                 communications for all passes
      */
      void receiveCommunication( const TotalArgumentType& totalArg ) const
      {
        // this method is called on the last pass which needs no
        // finalizing of communication, so call the correct method
        // on the previous pass
        // std::get<0> ( totalArg ) is the global argument type
        previousPass_.receiveCommunication( *(std::get<0> ( totalArg )) );
      }

      /** \brief initializeCommunication of this pass, this will initialize
                 the communication of destination_ and has to be overloaded in
                 the implementation
      */
      virtual void initComm() const {}

      /** \brief finalizeCommunication of this pass, this will collect
                 the communication of destination_ and has to be overloaded in
                 the implementation
      */
      virtual void finalizeComm() const {}

      /** \brief receiveCommunication of this pass,
                 which will reset changes the communication
                 did to the destination_ and has to be overloaded in
                 the implementation
      */
      virtual void receiveComm() const {}

    protected:
      //! destination (might be set from outside)
      DestinationType* destination_;

      //! object to delete destination_
      DeleteHandlerType* deleteHandler_;

      // previous pass
      PreviousPassType& previousPass_;

      // current calculation time
      double time_;

      // flag whether we are the last pass, i.e. we have to finalize the communication
      mutable bool finalizeCommunication_ ;
    }; // end class Pass

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_PASS_HH
