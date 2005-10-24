#ifndef DUNE_PASS_HH
#define DUNE_PASS_HH

//- System includes

//- Dune includes
#include <dune/common/tuples.hh>
#include <dune/common/operator.hh>

namespace Dune {

  /**
   * @brief End marker for a compile-time list of passes.
   *
   */
  template <class ArgumentImp>
  class StartPass {
  public:
    //- Enums and typedefs
    //! The start pass has index 0.
    enum {passNum=0};
    typedef ArgumentImp GlobalArgumentType;
    typedef Nil NextArgumentType;
    
  public:
    //- Public methods
    //! The pass method does nothing.
    void pass(const GlobalArgumentType& arg) const {}
    //! Returns the closure of the destination tuple.
    NextArgumentType localArgument() const { return Nil(); }
    //! No memory needs to be allocated.
    void allocateLocalMemory() {}
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
  template <class PassTraits, class PreviousPassImp>
  class Pass :
    public Operator<typename PreviousPassImp::GlobalArgumentType::DofType, 
                    typename PassTraits::DestinationType::DofType,
                    typename PreviousPassImp::GlobalArgumentType, 
                    typename PassTraits::DestinationType>
  {
    template <class PT, class PP>
    friend class Pass;

  public:
    //- Enums and typedefs
    //! The index of the pass.
    enum {passNum = PreviousPassImp::passNum + 1};

    //! Traits class provided by the implementor.
    typedef PassTraits Traits;
    //! Type of the preceding pass.
    typedef PreviousPassImp PreviousPassType;
    //! Type of the discrete function which stores the result of this pass' 
    //! computations.
    typedef typename PassTraits::DestinationType DestinationType;
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
    //- Public methods
    //! Constructor
    //! \param pass Previous pass
    Pass(PreviousPassType& pass) :
      destination_(0),
      previousPass_(pass)
    {
      // this ensures that the last pass doesn't allocate temporary memory
      // (its destination discrete function is provided from outside)
      previousPass_.allocateLocalMemory();
    }

    virtual ~Pass() {
      delete destination_;
      destination_ = 0;
    }

    //! \brief Application operator.
    //! The application operator is called by the client directly. It makes
    //! only sense to call this operator directly on the last pass.
    void operator()(const GlobalArgumentType& arg, DestinationType& dest) const
    {
      previousPass_.pass(arg);
      const TotalArgumentType totalArg(&arg, previousPass_.localArgument());
      this->compute(totalArg, dest);
    }

    //! Allocates the local memory of a pass, if needed.
    virtual void allocateLocalMemory() = 0;

  protected:
    // ? Really do this? Can't you allocate the memory directly?
    DestinationType* destination_;

  private:
    //! Same as application operator, but uses own memory instead of the
    //! discrete function provided by the client. This method is called on all
    //! passes except the last one.
    void pass(const GlobalArgumentType& arg) const {
      operator()(arg, *destination_);
    }

    //! Returns a compilation of the results of the preceding passes
    NextArgumentType localArgument() const {
      return NextArgumentType(destination_, previousPass_.localArgument());
    }

  private:
    //! Does the actual computations. Needs to be overridden in the derived
    //! clases
    virtual void compute(const TotalArgumentType& arg, 
                         DestinationType& dest) const = 0;

  private:
    PreviousPassType& previousPass_;
  }; // end class Pass


  //! Specialisation of Pass which provides a grid walk-through, but leaves
  //! open what needs to be done on each elements.
  template <class PassTraits, class PreviousPassImp>
  class LocalPass :
    public Pass<PassTraits, PreviousPassImp>
  {
  public:
    typedef PreviousPassImp PreviousPassType;

    typedef Pass<PassTraits, PreviousPassImp> BaseType;
    typedef typename BaseType::TotalArgumentType ArgumentType;

    typedef typename PassTraits::DestinationType DestinationType;
    typedef typename PassTraits::SpaceType SpaceType;
    typedef typename SpaceType::IteratorType IteratorType;
    typedef typename IteratorType::Entity Entity;

  public:
    //! Constructor
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function of this pass. 
    LocalPass(PreviousPassImp& pass, SpaceType& spc) :
      BaseType(pass),
      spc_(spc)
    {}

    virtual ~LocalPass() {}

    //! Build up local memory.
    virtual void allocateLocalMemory() 
    {
      if (!this->destination_) {
        this->destination_ = new DestinationType("temp", spc_);
      }
    }

  private:
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

  private:
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
    SpaceType& spc_;
  };
  
} // end namespace Dune

#endif
