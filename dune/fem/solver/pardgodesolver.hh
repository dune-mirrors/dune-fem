#ifndef PARDG_ODE_SOLVER_HH
#define PARDG_ODE_SOLVER_HH

#error DO NOT USE : odesolver.hh is the right one

// inlcude all used headers before, that they don not appear in DuneODE 

//- system includes 
#include <iostream>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cassert>
#include <sys/times.h>
#if HAVE_MPI
#include <mpi.h>
#endif

//- Dune includes 
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

//- include runge kutta ode solver 
#include <dune/fem/solver/rungekutta.hh>

// include headers of PARDG 
#include "pardg.hh"

namespace DuneODE {

#ifdef USE_PARDG_ODE_SOLVER

struct OdeSolverParameters
: public LocalParameter< OdeSolverParameters, OdeSolverParameters >
{
  virtual double nonLinearEps ( ) const
  {
    return Parameter::getValue< double >( "fem.ode.nonlineareps", 1e-5 );
  }
  virtual double linearEps ( ) const
  {
    return Parameter::getValue< double >( "fem.ode.lineareps", 1e-7 );
  }
  virtual bool nonLinearVerbose ( ) const
  {
    return (Parameter::getValue< int >( "fem.ode.nonlinearverbose", 0 ) == 1);
  }
  virtual bool linearVerbose ( ) const
  {
    return (Parameter::getValue< int >( "fem.ode.linearverbose", 0 ) == 1);
  }
  virtual int linearSolver ( ) const
  {
    static const std::string solvers[]
      = { "gmres", "bicgstab", "cg" };
    return Parameter::getEnum( "fem.ode.linearsolver", solvers, 0 );
  }
  virtual int gmresCycles ( ) const
  {
    return Parameter::getValue< int >( "fem.ode.gmresCycles", 15 );
  }
  virtual double cflSigma ( ) const
  {
    return Parameter::getValue< double >( "fem.ode.cflsigma", 1.1 );
  }
  virtual int cflUpperIterations ( ) const
  {
    return Parameter::getValue< int >( "fem.ode.cflupperiterations", 20 );
  }
  virtual int cflLowerIterations ( ) const
  {
    return Parameter::getValue< int >( "fem.ode.cflloweriterations", 10);
  }
  virtual bool cflFactor ( double iter,
                           double lowerIter, double upperIter,
                           double sigma,
                           double &cfl ) const
  {
    if ( iter < lowerIter ) 
    {
      cfl *= sigma;
      return true;
    }
    else if ( iter > upperIter ) 
    {
      cfl *= (double)upperIter/(sigma*(double)iter);
      return true;
    }
    return false;
  }
};


template <class Operator>
class OperatorWrapper : public PARDG::Function 
{
  // type of discrete function 
  typedef typename Operator::DestinationType DestinationType;
  // type of discrete function space 
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
 public:
  //! constructor 
  OperatorWrapper(Operator& op) 
    : op_(op) , space_(op_.space()) 
  {}

  //! apply operator applies space operator and creates temporary
  //! discrete function using the memory from outside 
  void operator()(const double *u, double *f, int i = 0) 
  {
    // create fake argument 
    DestinationType arg("ARG",space_,u);
    // create fake destination 
    DestinationType dest("DEST",space_,f);
    
    // set actual time of iteration step
    op_.setTime( this->time() );
    
    // call operator apply 
    op_(arg,dest);
  }

  //! return size of argument 
  int dim_of_argument(int i = 0) const 
  { 
    if (i==0) return op_.space().size();
    else 
    {
      assert(0);
      abort();
      return -1;
    }
  }
  
  //! return size of destination  
  int dim_of_value(int i = 0) const 
  { 
    if (i==0) return op_.space().size();
    else 
    {
      assert(0);
      abort();
      return -1;
    }
  }
private:
  // operator to call 
  Operator& op_;
  // discrete function space 
  const SpaceType& space_;
};


/**
   @ingroup ODESolver
   @{
 **/

/* \brief Explicit ODE Solver base class */
template<class Operator>
class ExplTimeStepperBase 
{
  typedef typename Operator::DestinationType DestinationType; 
public:
  ExplTimeStepperBase(Operator& op, 
                      Dune::TimeProviderBase& tp, 
                      int pord, 
                      bool verbose) :
    ord_(pord),
    comm_(PARDG::Communicator::instance()),
    op_(op),
    expl_(op),
    ode_(0),
    initialized_(false)
  {
    switch (pord) {
      case 1: ode_ = new PARDG::ExplicitEuler(comm_,expl_); break;
      case 2: ode_ = new PARDG::ExplicitTVD2(comm_,expl_); break;
      case 3: ode_ = new PARDG::ExplicitTVD3(comm_,expl_); break;
      case 4: ode_ = new PARDG::ExplicitRK4(comm_,expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }

    if (verbose)
    {
      ode_->DynamicalObject::set_output(cout);
    }
  }
 
  // initialize time step size 
  bool initialize (const DestinationType& U0) 
  {
    // initialized dt on first call
    if ( ! initialized_ )     
    {
      DestinationType tmp(U0);
      this->op_(U0,tmp);
      initialized_ = true;
      return true;
    }
    return false;
  }

  //! destructor  
  ~ExplTimeStepperBase() { delete ode_; }

  //! return reference to ode solver 
  PARDG::ODESolver& odeSolver() 
  {
    assert( ode_ );
    return *ode_;
  }
  void printmyInfo(string filename) const 
  {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "Explicit ODE solver, steps: " << this->ord_ << "\n\n";
    ofs.close();
    //this->op_.printmyInfo(filename);
  }
  
protected:
  int ord_;
  PARDG::Communicator & comm_;
  const Operator& op_;
  OperatorWrapper<Operator> expl_;
  PARDG::ODESolver* ode_;
  bool initialized_;
};

//! ExplicitOdeSolver 
template<class DestinationImp>
class ExplicitOdeSolver : 
  public OdeSolverInterface<DestinationImp> ,
  public ExplTimeStepperBase<SpaceOperatorInterface<DestinationImp> >  
{
  typedef DestinationImp DestinationType; 
  typedef SpaceOperatorInterface<DestinationType> OperatorType;
  typedef ExplTimeStepperBase<OperatorType> BaseType; 
public:
  //! constructor 
  ExplicitOdeSolver(OperatorType& op, Dune :: TimeProviderBase &tp, int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {}

  //! destructor 
  virtual ~ExplicitOdeSolver() {}
 
  //! initialize solver 
  void initialize(const DestinationType& U0)
  {
    BaseType :: initialize (U0);
    timeProvider_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
  }

  //! solve system 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! this->initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ExplicitOdeSolver wasn't initialized before first call!");
    }
    
    // get dt 
    const double dt = timeProvider_.deltaT();
    
    // should be larger then zero 
    assert( dt > 0.0 );
    
    // get time 
    const double time = timeProvider_.time();

    // get leakPointer 
    double* u = U0.leakPointer();
    
    // call ode solver 
    const bool convergence = this->odeSolver().step(time, dt , u);

    // set time step estimate of operator 
    timeProvider_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
    
    assert(convergence);
    if(!convergence) 
    {
      timeProvider_.invalidateTimeStep();
      std::cerr << "No Convergence of ExplicitOdeSolver! \n";
    }
  }

protected:
  Dune::TimeProviderBase& timeProvider_;
};

//////////////////////////////////////////////////////////////////
//
//  --ImplTimeStepperBase
//
//////////////////////////////////////////////////////////////////
template<class Operator>
class ImplTimeStepperBase
  : public OdeSolverInterface<typename Operator::DestinationType> 
{
  typedef typename Operator :: DestinationType DestinationType; 
public:
  template <class ParDGODESolver>
  ImplTimeStepperBase(const Operator &op,
                      ParDGODESolver* ode,
                      Dune :: TimeProviderBase &tp, 
                      const OdeSolverParameters& parameter=OdeSolverParameters() )
  : param_(parameter.clone()),
    op_(op),
    ode_(ode),
    timeProvider_(tp),
    comm_(PARDG::Communicator::instance()),
    linsolver_(0),
    minIter_( param_->cflLowerIterations() ),
    maxIter_( param_->cflUpperIterations() ),
    sigma_( param_->cflSigma() ),
    cfl_(1),
    initialized_(false)
  {
    switch (param_->linearSolver())
    {
      case 0: linsolver_ = new PARDG::GMRES( comm_,param_->gmresCycles() ); break;
      case 1: linsolver_ = new PARDG::BICGSTAB( comm_ ); break;
      case 2: linsolver_ = new PARDG::CG( comm_ ); break;
    };
    linsolver_->set_tolerance(param_->linearEps(),false);
    linsolver_->set_max_number_of_iterations(10000);
    ode->set_linear_solver(*linsolver_);
    ode->set_tolerance(param_->nonLinearEps());
    ode->set_max_number_of_iterations(150);
    
    if( param_->linearVerbose() ) 
    {
      linsolver_->IterativeSolver::set_output(cout);
      linsolver_->IterativeSolver::set_output(cout);
    }
    if( param_->nonLinearVerbose() ) 
    {
      ode_->DynamicalObject::set_output(cout);
      ode_->DynamicalObject::set_output(cout);
    }
  }

  //! destructor 
  ~ImplTimeStepperBase() {delete ode_; delete linsolver_;}
  
  // initialize time step size 
  void initialize (const DestinationType& U0) 
  {
    // initialized dt on first call
    if ( ! initialized_ )     
    {
      DestinationType tmp(U0);
      op_(U0,tmp);
      initialized_ = true;
      // set time step estimate of operator 
      timeProvider_.provideTimeStepEstimate( cfl_ * op_.timeStepEstimate() );
    }
  }

  //! solve 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"Implicit Ode solver wasn't initialized before first call!");
    }

    const double dt   = timeProvider_.deltaT();
    assert( dt > 0.0 );
    const double time = timeProvider_.time();

    // get pointer to solution
    double* u = U0.leakPointer();
      
    const bool convergence = ode_->step(time , dt , u);
    const int iter = linsolver_->number_of_iterations();
    // set time step estimate of operator 
    if (convergence) 
    {
      // control the number of iterations of the linear solver
      // the values for min_it and max_it has to be determined by experience
      if ( param_->cflFactor( iter, minIter_, maxIter_,sigma_, cfl_) 
           && Parameter::verbose() )
      {
        derr << "ODESolver: New cfl number is: "<< cfl_ << "\n";
        derr << "           number of iterations of linear solver  " << iter << std::endl;
      }
      
      timeProvider_.provideTimeStepEstimate( cfl_ * op_.timeStepEstimate() );
    
      linsolver_->reset_number_of_iterations();
    }
    else 
    {
      cfl_ *= 0.5;
      timeProvider_.provideTimeStepEstimate( cfl_ * dt );
      timeProvider_.invalidateTimeStep();
     
      // output only in verbose mode 
      if( Parameter :: verbose () )
      {
        derr << "ODE solver did not convergence - new cfl number is "<< cfl_ << std :: endl;
      }
    }
  }
  void printmyInfo(string filename) const {
  }
  
protected:
  const OdeSolverParameters* param_;
  const Operator& op_;
  PARDG::ODESolver* ode_;
  Dune :: TimeProviderBase &timeProvider_;
  PARDG::Communicator & comm_;   
  PARDG::IterativeLinearSolver* linsolver_;
  double minIter_,maxIter_,sigma_;
  double cfl_;
  bool initialized_;
};


///////////////////////////////////////////////////////
//
//  --ImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class ImplicitOdeSolver : 
  public ImplTimeStepperBase<SpaceOperatorInterface<DestinationImp> > 
{
  typedef DestinationImp DestinationType;
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef ImplTimeStepperBase<OperatorType> BaseType;

public:
  ImplicitOdeSolver(OperatorType& op, 
                    Dune::TimeProviderBase& tp,
                    int pord, 
                    const OdeSolverParameters& parameter=OdeSolverParameters() )
  : BaseType(op, setOdeSolver(pord,PARDG::Communicator::instance()),tp,parameter),
    impl_(op)
  { 
  }

  virtual ~ImplicitOdeSolver() {}

private:
  PARDG::DIRK *setOdeSolver(int order, PARDG::Communicator& comm)
  {
    switch (order) 
    {
      case 1: return new PARDG::ImplicitEuler(comm,impl_);
      case 2: return new PARDG::Gauss2(comm,impl_);
      case 3: return new PARDG::DIRK3(comm,impl_);
      //case 4: ode_ = new PARDG::ExplicitRK4(comm,expl_); break;
      default : std::cerr << "DI-Runge-Kutta method of this order not implemented in PARDG package" 
                          << std::endl;
                abort();
    }
  }
  OperatorWrapper<OperatorType> impl_;
}; // end ImplicitOdeSolver

///////////////////////////////////////////////////////
//
//  --SemiImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class SemiImplicitOdeSolver : 
  public ImplTimeStepperBase<SpaceOperatorInterface<DestinationImp> > 
{
  typedef DestinationImp DestinationType;
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef ImplTimeStepperBase<OperatorType> BaseType;

public:
  SemiImplicitOdeSolver(OperatorType& opExpl,
                        OperatorType& opImpl,
                        Dune::TimeProviderBase& tp,
                        int pord, 
                        const OdeSolverParameters& parameter=OdeSolverParameters() )
  : BaseType(opExpl, setOdeSolver(pord,PARDG::Communicator::instance()),tp,parameter),
    expl_(opExpl),
    impl_(opImpl)
  { 
  }

  virtual ~SemiImplicitOdeSolver() {}

private:
  PARDG::SIRK *setOdeSolver(int order,PARDG::Communicator& comm)
  {
    switch (order) 
    {
      case 1: return new PARDG::SemiImplicitEuler(comm,impl_,expl_); 
      case 2: return new PARDG::IMEX_SSP222(comm,impl_,expl_); 
      case 3: return new PARDG::SIRK33(comm,impl_,expl_); 
      default : std::cerr << "IMEX-Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }
  }
  OperatorWrapper<OperatorType> expl_;
  OperatorWrapper<OperatorType> impl_;
}; // end SemiImplicitOdeSolver

///////////////////////////////////////////////////////
//
//  --general OdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class OdeSolver : 
  public OdeSolverInterface<DestinationImp>
{
  typedef DestinationImp DestinationType;
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;

public:
  OdeSolver(OperatorType& explOp, 
            OperatorType& implOp, 
            Dune::TimeProviderBase& tp,
            int pord, 
            const OdeSolverParameters& parameter=OdeSolverParameters() )
  {
    create( &explOp, &implOp, tp, pord, parameter );
  }
  OdeSolver(OperatorType& explOp, 
            OperatorType* implOp, 
            Dune::TimeProviderBase& tp,
            int pord, 
            const OdeSolverParameters& parameter=OdeSolverParameters() )
  {
    create( &explOp, implOp, tp, pord, parameter );
  }
  OdeSolver(OperatorType* explOp, 
            OperatorType& implOp, 
            Dune::TimeProviderBase& tp,
            int pord, 
            const OdeSolverParameters& parameter=OdeSolverParameters() )
  {
    create( explOp, &implOp, tp, pord, parameter );
  }
  OdeSolver(OperatorType* explOp, 
            OperatorType* implOp, 
            Dune::TimeProviderBase& tp,
            int pord, 
            const OdeSolverParameters& parameter=OdeSolverParameters() )
  {
    create( explOp, implOp, tp, pord, parameter );
  }

  virtual ~OdeSolver() 
  {
    delete ode_;
  }
  
  //! initialize solver 
  void initialize(const DestinationType& U0)
  {
    // initialize solver 
    ode_->initialize (U0);
  }

  //! solve 
  void solve(DestinationType& U0) 
  {
    ode_->solve(U0);
  }
private:
  void create(OperatorType* explOp, 
              OperatorType* implOp, 
              Dune::TimeProviderBase& tp,
              int pord, 
              const OdeSolverParameters& parameter )
  {
    assert( explOp || implOp );
    if (implOp == 0)
      ode_ = new ExplicitOdeSolver<DestinationType>(*explOp,tp,pord,false);
    else if (explOp == 0)
      ode_ = new ImplicitOdeSolver<DestinationType>(*implOp,tp,pord,parameter);
    else
      ode_ = new SemiImplicitOdeSolver<DestinationType>(*explOp,*implOp,tp,pord,parameter);
  }
  OdeSolverInterface<DestinationImp> *ode_;
};

#endif // USE_PARDG_ODE_SOLVER
 
/**
 @} 
**/

} // end namespace DuneODE

#undef USE_EXTERNAL_BLAS
#endif
