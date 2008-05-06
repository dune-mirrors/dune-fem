#ifndef DUNE_PASS_HH
#define DUNE_PASS_HH

//- System includes
#include <string>
#include <limits>

//- Dune includes
#include <dune/common/timer.hh>
#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh>

//- local includes 
#include <dune/fem/operator/common/operator.hh>

namespace Dune {

  /*! @addtogroup Pass 
   *
   */
  /*!
   * @brief End marker for a compile-time list of passes.
   *
   */
  template < class ArgumentImp , int passIdImp  = -1 >
  class StartPass {
  public:
    //- Enums and typedefs
    //! The start pass has index 0.
    enum {passNum=0};
    enum{ passId = passIdImp };
    //! The argument (and destination) type of the overall operator
    typedef ArgumentImp GlobalArgumentType;
    //! End marker for tuple of return types of the passes
    typedef Nil NextArgumentType;
    
  public:
    //! empty constructor 
    StartPass() {} 
    //! copy constructor 
    StartPass(const StartPass&) {}
    
    //- Public methods
    //! The pass method does nothing.
    void pass(const GlobalArgumentType& arg) const {
    }
    void printTexInfo(std::ostream& out) const {
    }
    //! Returns the closure of the destination tuple.
    NextArgumentType localArgument() const { return nullType(); }
    //! No memory needs to be allocated.
    void allocateLocalMemory() {}

    //! StartPass does not need a time 
    void setTime(const double) {}

    //! time step estimate here is max 
    double timeStepEstimateImpl () const 
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
    public Operator<typename PreviousPassImp::GlobalArgumentType::RangeFieldType, 
                    typename DiscreteModelImp::Traits::DestinationType::RangeFieldType,
                    typename PreviousPassImp::GlobalArgumentType, 
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

    //! Type of the discrete function which is passed to the overall operator
    //! by the user
    typedef typename PreviousPassType::GlobalArgumentType GlobalArgumentType;
    //! Tuple containing destination types of all preceding passes.
    typedef typename PreviousPassType::NextArgumentType LocalArgumentType;
    //! Tuple containing destination types of all preceding passes plus the
    //! global argument. This serves as the argument for this pass' 
    //! computations
    typedef Pair<
      const GlobalArgumentType*, LocalArgumentType> TotalArgumentType;
    //! Tuple containing destination types of all passes up to this one.
    typedef Pair<DestinationType*, LocalArgumentType> NextArgumentType;

  public:
    // * Hack!!!! Remove again
    int passNumber() const { return passNum; }

    //- Public methods
    //! Constructor
    //! \param pass Previous pass
    Pass(PreviousPassType& pass) :
      destination_(0),
      deleteHandler_(0),
      previousPass_(pass),
      time_(0.0)
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
      typename PreviousPassType::NextArgumentType prevArg=previousPass_.localArgument();
      const TotalArgumentType totalArg(&arg, prevArg);
      this->compute(totalArg, dest);
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
      return std::min( previousPass_.timeStepEstimateImpl(),
                       this->timeStepEstimateImpl() );
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
    void pass(const GlobalArgumentType& arg) const {
      operator()(arg, *destination_);
    }

    //! Returns a compilation of the results of the preceding passes
    NextArgumentType localArgument() const {
      typename PreviousPassType :: NextArgumentType nextArg( previousPass_.localArgument() );
      return NextArgumentType(destination_, nextArg );
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

  protected:
    // ? Really do this? Can't you allocate the memory directly?
    DestinationType* destination_;
    //! object to delete destination_ 
    DeleteHandlerType* deleteHandler_; 

    // previous pass 
    PreviousPassType& previousPass_;

    // current calculation time 
    double time_;
  }; // end class Pass


  //! Specialisation of Pass which provides a grid walk-through, but leaves
  //! open what needs to be done on each elements.
  template <class DiscreteModelImp, class PreviousPassImp , int passIdImp >
  class LocalPass :
    public Pass< DiscreteModelImp , PreviousPassImp , passIdImp>
    //public Pass<DiscreteModelImp, PreviousPassImp>::DeleteHandler
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
    typedef typename IteratorType::Entity Entity;

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
    virtual void applyLocal(Entity& en) const = 0;

  protected:
    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const 
    {
      // get stopwatch 
      Timer timer; 
      
      prepare(arg, dest);
      
      IteratorType endit = spc_.end();
      for (IteratorType it = spc_.begin(); it != endit; ++it) {
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
  
} // end namespace Dune

#endif
