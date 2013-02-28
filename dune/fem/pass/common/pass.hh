#ifndef DUNE_FEM_PASS_HH
#define DUNE_FEM_PASS_HH

#include <limits>
#include <string>

#include <dune/common/timer.hh>
#include <dune/common/tuples.hh>

#include <dune/fem/operator/common/operator.hh>

#include "tupleutility.hh"

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
    template < class ArgumentImp , int passIdImp  = -1,
               class NonBlockingCommunication = EmptyNonBlockingComm  >
    class StartPass 
    {
      NonBlockingCommunication nonBlockingComm_;
    public:
      //- Enums and typedefs
      //! The start pass has index 0.
      enum{ passNum=0 };
      enum{ passId = passIdImp };
      //! The argument (and destination) type of the overall operator
      typedef ArgumentImp GlobalArgumentType;
      //! End marker for tuple of return types of the passes
      typedef Dune::tuple<> NextArgumentType;
      
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

      // dummy printTexInfo method 
      void printTexInfo(std::ostream& out) const {}

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
      enum{ passId = passIdImp };
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
        static DeleteHandler<ObjectToDelete>& instance () 
        {
          static DeleteHandler<ObjectToDelete> mh;
          return mh;
        }
      };

      //- Enums and typedefs
      //! The index of the pass.
      enum {passNum = PreviousPassImp::passNum + 1};

      //! Type of the preceding pass.
      typedef PreviousPassImp PreviousPassType;
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
        typename PreviousPassType::NextArgumentType prevArg = previousPass_.localArgument();
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
      NextArgumentType localArgument() const {
        typename PreviousPassType :: NextArgumentType nextArg( previousPass_.localArgument() );
        return NextArgumentType(destination_, nextArg );
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

    protected:
      //! Does the actual computations. Needs to be overridden in the derived
      //! clases
      virtual void compute(const TotalArgumentType& arg, 
                           DestinationType& dest) const = 0;

      //! derived passes have to implement this method 
      //! returning the time step estimate 
      virtual double timeStepEstimateImpl() const 
      {
        return std::numeric_limits<double>::max(); 
      }

      /** \brief finalizeCommunication collects possbily initiated non-blocking
                 communications for all passes 
      */
      void finalizeCommunication( const TotalArgumentType& totalArg ) const 
      {
        // this method is called on the last pass which needs no 
        // finalizing of communication, so call the correct method 
        // on the previous pass 
        // totalArg.first() is the global argument type 
        previousPass_.finalizeCommunication( *totalArg.first() );
      }

      /** \brief receiveCommunication collects possbily initiated non-blocking
                 communications for all passes 
      */
      void receiveCommunication( const TotalArgumentType& totalArg ) const 
      {
        // this method is called on the last pass which needs no 
        // finalizing of communication, so call the correct method 
        // on the previous pass 
        // totalArg.first() is the global argument type 
        previousPass_.receiveCommunication( *totalArg.first() );
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


    //! Specialisation of Pass which provides a grid walk-through, but leaves
    //! open what needs to be done on each elements.
    template <class DiscreteModelImp, class PreviousPassImp , int passIdImp >
    class LocalPass :
      public Pass< DiscreteModelImp , PreviousPassImp , passIdImp>
    {
    public:
      //! Type of the preceding pass
      typedef PreviousPassImp PreviousPassType;

      //! Base class
      typedef Pass< DiscreteModelImp , PreviousPassImp , passIdImp > BaseType;
      //! The type of the argument (and destination) type of the overall
      //! operator
      
      typedef typename BaseType::TotalArgumentType ArgumentType;

      //! The discrete function representing the return value of this pass
      typedef typename DiscreteModelImp::Traits::DestinationType DestinationType;
      //! The discrete function space belonging to DestinationType
      typedef typename DiscreteModelImp::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      //! Iterator over the space
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      //! The codim 0 entity
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType ;

      // deprecated type 
      typedef EntityType  Entity;

    public:
      //! Constructor
      //! \param pass Previous pass
      //! \param spc Space belonging to the discrete function of this pass. 
      //! \param passName an identifier for this pass
      LocalPass(PreviousPassImp& pass, const DiscreteFunctionSpaceType& spc,
                std::string passName = "LocalPass") :
        BaseType(pass),
        spc_(spc),
        passName_(passName ),
        computeTime_(0.0)
      {}

      //! Destructor
      virtual ~LocalPass() {}

      //! Build up local memory.
      virtual void allocateLocalMemory() 
      {
        if (!this->destination_) 
        {
          std::ostringstream funcName;
          funcName << passName_ << "_" << this->passNumber();
          this->destination_ = new DestinationType(funcName.str(), spc_);
          
          // set mem handle for deleting destination_ 
          this->deleteHandler_ = &(BaseType::DeleteHandlerType::instance());
        }
      }

      //! return reference to space 
      const DiscreteFunctionSpaceType& space() const 
      { 
        return spc_; 
      }

      //! return accumulated time needed by pass's operator ()  
      //! this method also resets the compute time to zero 
      virtual double computeTime() const 
      {
        double ct = computeTime_;
        computeTime_ = 0.0;
        return ct;
      }

    protected:
      //! Actions to be carried out before a global grid walkthrough.
      //! To be overridden in a derived class.
      virtual void prepare(const ArgumentType& arg, 
                           DestinationType& dest) const = 0;
      //! Actions to be carried out after a global grid walkthrough. 
      //! To be overridden in a derived class.
      virtual void finalize(const ArgumentType& arg, 
                            DestinationType& dest) const = 0;
      //! Actions to be taken on every element. To be overridden in a derived 
      //! class.
      virtual void applyLocal( const EntityType& en ) const = 0;

    protected:
      //! The actual computations are performed as follows. First, prepare
      //! the grid walkthrough, then call applyLocal on each entity and then
      //! call finalize.
      void compute(const ArgumentType& arg, DestinationType& dest) const 
      {
        // get stopwatch 
        Dune::Timer timer; 
        
        prepare(arg, dest);
        
        IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
          applyLocal(*it);
        }

        finalize(arg, dest);

        // accumulate time 
        computeTime_ += timer.elapsed();
      }    

    protected:
      const DiscreteFunctionSpaceType& spc_;
      const std::string passName_;
      mutable double computeTime_;
    };
    
  } // namespace Fem 

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: Pass ;
  using Fem :: LocalPass ;
  using Fem :: StartPass ;

#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_HH
