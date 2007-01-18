#ifndef DUNE_PASS_HH
#define DUNE_PASS_HH

//- System includes
#include <string>

//- Dune includes
#include <dune/common/tuples.hh>

//- local includes 
#include <dune/fem/operator/common/operator.hh>

// * must go away
#include <dune/fem/misc/timeutility.hh>

namespace Dune {

  /*! @defgroup Pass Pass concept
   *  @ingroup OperatorCommon
   *
   */
  /*!
   * @brief End marker for a compile-time list of passes.
   *
   */
  template <class ArgumentImp>
  class StartPass {
  public:
    //- Enums and typedefs
    //! The start pass has index 0.
    enum {passNum=0};
    //! The argument (and destination) type of the overall operator
    typedef ArgumentImp GlobalArgumentType;
    //! End marker for tuple of return types of the passes
    typedef Nil NextArgumentType;
    
  public:
    //- Public methods
    //! The pass method does nothing.
    void pass(const GlobalArgumentType& arg) const {
    }
    //! Returns the closure of the destination tuple.
    NextArgumentType localArgument() const { return nullType(); }
    //! No memory needs to be allocated.
    void allocateLocalMemory() {}
    //! We don't need no time either
    void timeProvider(TimeProvider*) {}
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
  template <class DiscreteModelImp, class PreviousPassImp>
  class Pass :
    public Operator<typename PreviousPassImp::GlobalArgumentType::DofType, 
                    typename DiscreteModelImp::Traits::DestinationType::DofType,
                    typename PreviousPassImp::GlobalArgumentType, 
                    typename DiscreteModelImp::Traits::DestinationType>
  {
    template <class PT, class PP>
    friend class Pass;
  public:
    //! little interface class for deleting discrete function 
    //! held by this class 
    template <class ObjectToDelete> 
    struct MemHandler
    {
      virtual void freeLocalMemory(ObjectToDelete * obj) = 0;
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
    typedef MemHandler<DestinationType> MemHandlerType;

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
      memHandler_(0),
      previousPass_(pass)
    {
      // this ensures that the last pass doesn't allocate temporary memory
      // (its destination discrete function is provided from outside).
      previousPass_.allocateLocalMemory();
    }

    //! Destructor
    virtual ~Pass() 
    {
      // if memHandler was set by derived class, then use to delete destination_ 
      if( memHandler_ ) memHandler_->freeLocalMemory(destination_);
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
    virtual void allocateLocalMemory() = 0;

    //! Set time provider (which gives you access to the global time).
    void timeProvider(TimeProvider* time) {
      previousPass_.timeProvider(time);
      processTimeProvider(time);
    }

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

    //! With this method, you can get access to the TimeProvider, if needed.
    //! By default, nothing is done, ie the TimeProvider is discarded.
    virtual void processTimeProvider(TimeProvider* time) {}

  protected:
    //! Does the actual computations. Needs to be overridden in the derived
    //! clases
    virtual void compute(const TotalArgumentType& arg, 
                         DestinationType& dest) const = 0;

  protected:
    // ? Really do this? Can't you allocate the memory directly?
    DestinationType* destination_;
    //! object to delete destination_ 
    MemHandlerType * memHandler_; 

    // previous pass 
    PreviousPassType& previousPass_;
    
  }; // end class Pass


  //! Specialisation of Pass which provides a grid walk-through, but leaves
  //! open what needs to be done on each elements.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalPass :
    public Pass<DiscreteModelImp, PreviousPassImp>,
    public Pass<DiscreteModelImp, PreviousPassImp>::template 
              MemHandler<typename DiscreteModelImp::Traits::DestinationType>
  {
  public:
    //! Type of the preceding pass
    typedef PreviousPassImp PreviousPassType;

    //! Base class
    typedef Pass<DiscreteModelImp, PreviousPassImp> BaseType;
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
      passName_(passName )
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
        // set memHandler for deleting destination_ 
        this->memHandler_ = this;
      }
    }

    //! deletes destination_ 
    virtual void freeLocalMemory(DestinationType * obj) 
    {
      delete this->destination_;
      this->destination_ = 0;
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
      prepare(arg, dest);
      
      IteratorType endit = spc_.end();
      for (IteratorType it = spc_.begin(); it != endit; ++it) {
        applyLocal(*it);
      }

      finalize(arg, dest);
    }    

  private:
    const DiscreteFunctionSpaceType& spc_;
    const std::string passName_;
  };
  
} // end namespace Dune

#endif
