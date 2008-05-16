#ifndef DUNE_ELLIPTPASS_HH
#define DUNE_ELLIPTPASS_HH

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

//- local includes 
#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/ellipticdiscretemodel.hh>
#include <dune/fem/pass/ellipticmodelcaller.hh>

#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune {
/*! @addtogroup PassEllipt
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/

  //! Concrete implementation of Pass for elliptic equations using LDG
  template< class DiscreteModelImp , class PreviousPassImp , int passId = -1 >
  class LocalDGElliptGradientPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp , passId > 
  {
    typedef LocalDGElliptGradientPass< DiscreteModelImp , PreviousPassImp , passId > ThisType;
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassImp , passId > BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef EllipticDiscreteModelCaller
      < DiscreteModelType, ArgumentType, CombinedSelectorType >
      DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType::polynomialOrder };

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LocalDGElliptGradientPass(DiscreteModelType& problem, 
                    PreviousPassType& pass, 
                    const DiscreteFunctionSpaceType& spc) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      spc_(spc)
    {
    }
    void printTexInfo(std::ostream& out) const {
      BaseType::printTexInfo(out);
      out << "LocalDGElliptGradientPass: ";
      out << "\\\\ \n";
    }

    //! don't allocate memory here 
    virtual void allocateLocalMemory() {}
   
    //! return reference to caller 
    DiscreteModelCallerType & caller () { return caller_; }

    //! return previous pass of this pass 
    PreviousPassType & previousPass() { return this->previousPass_; }

    //! return problem for real fe pass 
    DiscreteModelType & problem () { return problem_; }

    //! return reference to space 
    const DiscreteFunctionSpaceType & space () const { return spc_; }
   
    //! Destructor
    virtual ~LocalDGElliptGradientPass() {
    }

    //! do nothing here 
    void applyLocal(EntityType& en) const
    {
    }

    //! do nothing here  
    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      abort();
    }

    //! do nothing here 
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! Some timestep size management.
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    const DiscreteFunctionSpaceType& spc_;
  };

  
  //! Concrete implementation of Pass for LDG.
  template< class DiscreteModelImp , class PreviousPassImp , int passId = -1 >
  class LocalDGElliptPass :
    public LocalPass< DiscreteModelImp , typename PreviousPassImp:: PreviousPassType , passId > 
  {
    typedef LocalDGElliptPass< DiscreteModelImp , PreviousPassImp , passId > ThisType;
  public:
    typedef typename PreviousPassImp::PreviousPassType  PreviousPassType;
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassType , passId > BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp GradFePassImp;

    typedef typename GradFePassImp :: DestinationType GradDestinationType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType;
    
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType, SelectorType > CombinedSelectorType;
    typedef DiscreteModelCaller
      < DiscreteModelType, ArgumentType, CombinedSelectorType >
      DiscreteModelCallerType;
   
    // Range of the destination
    enum { dimRange  = DiscreteFunctionSpaceType :: dimRange }; 
    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain }; 
                    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };

    // define previous pass of grad pass as previous pass of hole ellipt
    // pass
    typedef typename GradFePassImp :: PreviousPassType ElliptPrevPassType;

    typedef typename DiscreteModelType :: Traits :: template 
      LocalOperatorSelector<GradFePassImp> LocalOperatorSelectorType;
    typedef typename  LocalOperatorSelectorType :: LocalOperatorType LocalOperatorType;
    typedef typename  LocalOperatorSelectorType :: InverseOperatorType InverseOperatorType;

    //! type of restrict and prolong operator during adaptation 
    typedef LocalOperatorType RestrictProlongOperatorType;

  private:  
    DiscreteModelType& problem_; 
    const DiscreteFunctionSpaceType& spc_;
    const bool verbose_;
    mutable LocalOperatorType op_;

    const double eps_;
    const int maxIterFactor_; 
    mutable int maxIter_;

    InverseOperatorType invOp_; 

    mutable DestinationType rhs_;

    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
    mutable CommunicationManagerType comm_;
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param paramFile file name of parameter file to read various variables 
    //! 
    //!  NOTE: parameter read by this class 
    //!         - InvSolverEps epsilon for interative solver, default is 1e-10 
    //!         - verbose if true some output is given, default is false
    LocalDGElliptPass(DiscreteModelType& problem, 
                PreviousPassImp & pass, 
                const DiscreteFunctionSpaceType& spc,
                const std::string paramFile = "")
      : BaseType(pass.previousPass(),spc)
      , problem_(problem)
      , spc_(spc) 
      , verbose_(readVerbose(paramFile, spc_.grid().comm().rank() == 0))
      , op_(problem,pass,pass.previousPass(),spc,paramFile)
      , eps_(readEps(paramFile, verbose_ ))
      , maxIterFactor_(4) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps_,eps_,maxIter_,verbose_)
      , rhs_("FEPass::RHS",spc)
      , comm_(spc_)
    {
      //assert( this->destination_ );
    }

     //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param dest Pointer to needed temporary memory (otherwise created by pass)
    //! \param paramFile file name of parameter file to read various variables 
    //! 
    //!  NOTE: parameter read by this class 
    //!         - InvSolverEps epsilon for interative solver, default is 1e-10 
    //!         - verbose if true some output is given, default is false  
    LocalDGElliptPass(DiscreteModelType& problem, 
                PreviousPassImp & pass, 
                DestinationType & dest,
                const std::string paramFile = "")
      : BaseType(pass.previousPass(),dest.space())
      , problem_(problem)
      , spc_(dest.space()) 
      , verbose_(readVerbose(paramFile, spc_.grid().comm().rank() == 0))
      , op_(problem,pass,pass.previousPass(),spc_,paramFile)
      , eps_(readEps(paramFile, verbose_ ))
      , maxIterFactor_(4) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps_,eps_,maxIter_,verbose_)
      , rhs_("FEPass::RHS",spc_)
      , comm_(spc_)
    {
      assert( this->destination_ == 0 );
      this->destination_ = &dest;
    }
    void printTexInfo(std::ostream& out) const {
      BaseType::printTexInfo(out);
      out << "LocalDGElliptPass: ";
      out << " eps = " << eps_
          << " inverse Operator: "
          << "\\\\ \n";
      // invOp_.printTexInfo(out);
    }

    //! do nothing here 
    void applyLocal(EntityType& en) const
    {
    }
    
    GradDestinationType & tmpMemory() { return op_.tmpMemory(); }

    //! return restrict and prolong operator for fe-pass 
    RestrictProlongOperatorType & restrictProlongOperator () { return op_; }

    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // re-compute matrix 
      op_.computeMatrix( arg, rhs_ );
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! compute method 
    virtual void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // prepare operator 
      prepare(arg,dest);

      // calculate new maxIter  
      maxIter_ = maxIterFactor_ * spc_.size();

      // solve the system 
      invOp_(rhs_,dest);

      // do data exchange 
      comm_.exchange( dest );
    } 

    template <class FuncType, class GradType>
    void evalGradient(const FuncType & u, GradType & grad, bool applyMass = false) const
    {
      op_.evalGradient(u,grad,applyMass);
    }
  private:
    bool readVerbose(const std::string& paramFile, const bool verboseOutput) const 
    {
      if( paramFile == "" ) return false;
      int val = 0;
      readParameter(paramFile,"verbose",val, verboseOutput);
      return ( (val == 1) && verboseOutput) ? true : false;
    }
    
    double readEps(const std::string& paramFile, const bool output) const 
    {
      double eps = 1e-10; 
      if( paramFile == "" ) return eps;
      readParameter(paramFile,"InvSolverEps",eps, output);
      return eps;
    }
  };

  //! Concrete implementation of Pass for LDG.
  template< class DiscreteModelImp , class PreviousPassImp , int passId = -1 >
  class LocalDGElliptGradPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp, passId > 
  {
    typedef LocalDGElliptGradPass< DiscreteModelImp , PreviousPassImp , passId > ThisType;
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassImp , passId > BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef DiscreteModelCaller
      < DiscreteModelType, ArgumentType, CombinedSelectorType >
      DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd =DiscreteFunctionSpaceType::polOrd};

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param factor Sign of the gradient
    //! \param applyMassMatrix if true also mass matrix (if exsists) 
    //! is applied when evaluating the gradient, default is true 
    LocalDGElliptGradPass(DiscreteModelType& problem, 
                    PreviousPassType& pass, 
                    const DiscreteFunctionSpaceType& spc,
                    double factor = -1.0,
                    bool applyMassMatrix = true) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      spc_(spc),
      prevPass_(pass),
      comm_(spc_),
      factor_(factor),
      applyMassMatrix_(applyMassMatrix)
    {
    }
    void printTexInfo(std::ostream& out) const {
      BaseType::printTexInfo(out);
      out << "LocalDGElliptGradPass: ";
      out << " factor = " << factor_
          << " apply Mass = " << applyMassMatrix_
          << "\\\\ \n";
    }

    //! don't allocate memory here 
    virtual void allocateLocalMemory() {}
   
    //! return reference to caller 
    DiscreteModelCallerType & caller () { return caller_; }

    //! return previous pass of this pass 
    PreviousPassType & previousPass() { return this->previousPass_; }

    //! return problem for real fe pass 
    DiscreteModelType & problem () { return problem_; }

    const DiscreteFunctionSpaceType & space () const { return spc_; }
   
    //! Destructor
    virtual ~LocalDGElliptGradPass() {
    }

    //! calls evalGradient of previous pass 
    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      // normal call procedure 
      prevPass_.pass(arg);

      // now get gradient from previous pass 
      prevPass_.evalGradient(prevPass_.destination(), dest,
                             applyMassMatrix_);

      // return -grad p 
      dest *= factor_;

      // exchange data 
      comm_.exchange( dest );
    }

    //! nothing to do here
    void applyLocal(EntityType& en) const
    {
    }
    
    //! nothing to do here
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! nothing to do here 
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    const DiscreteFunctionSpaceType& spc_;
    mutable PreviousPassImp & prevPass_;
    mutable CommunicationManagerType comm_;
    const double factor_;
    const bool applyMassMatrix_;
  };

} // end namespace Dune

#endif
