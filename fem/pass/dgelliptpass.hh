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

namespace Dune {
/*! @addtogroup PassEllipt
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/
  //! Concrete implementation of Pass for DG.
  template< class DiscreteModelImp , class PreviousPassImp , int passId = -1 >
  class LocalDGElliptPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp , passId > 
  {
    typedef LocalDGElliptPass< DiscreteModelImp , PreviousPassImp , passId > ThisType;
  public:
    typedef PreviousPassImp PreviousPassType;
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassType , passId > BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

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
    typedef DiscreteModelCallerDefault
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

    typedef typename DiscreteModelType :: Traits :: template 
      LocalOperatorSelector<PreviousPassType> LocalOperatorSelectorType;

    typedef typename  LocalOperatorSelectorType :: LocalOperatorType LocalOperatorType;
    typedef typename  LocalOperatorSelectorType :: InverseOperatorType InverseOperatorType;

    //! type of restrict and prolong operator during adaptation 
    typedef LocalOperatorType RestrictProlongOperatorType;

  protected:  
    DiscreteModelType& problem_; 
    const DiscreteFunctionSpaceType& spc_;
    const bool verbose_;
    mutable LocalOperatorType op_;

    const double eps_;
    const int maxIterFactor_; 
    mutable int maxIter_;

    InverseOperatorType invOp_; 

    mutable DestinationType rhs_;

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
      : BaseType(pass,spc)
      , problem_(problem)
      , spc_(spc) 
      , verbose_(readVerbose(paramFile, spc_.grid().comm().rank() == 0))
      , op_(problem,pass,spc,paramFile)
      , eps_(readEps(paramFile, verbose_ ))
      , maxIterFactor_(4) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps_,eps_,maxIter_,verbose_)
      , rhs_("DGPass::RHS",spc)
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
      : BaseType(pass,dest.space())
      , problem_(problem)
      , spc_(dest.space()) 
      , verbose_(readVerbose(paramFile, spc_.grid().comm().rank() == 0))
      , op_(problem,pass,spc_,paramFile)
      , eps_(readEps(paramFile, verbose_ ))
      , maxIterFactor_(4) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps_,eps_,maxIter_,verbose_)
      , rhs_("DGPass::RHS",spc_)
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

    //! Set time provider (which gives you access to the global time).
    void setTime(const double t)
    {
      BaseType :: setTime( t );
      op_.setTime( t );
    }

    //! do nothing here 
    void applyLocal( const EntityType& en) const
    {
    }
    
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
      spc_.communicate( dest );
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

} // end namespace Dune
#endif
