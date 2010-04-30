#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

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
struct ODEParameters
: public LocalParameter< ODEParameters, ODEParameters >
{ 
  ODEParameters() : 
    min_it( Parameter::getValue< int >( "fem.ode.miniterations" , 14 ) ),
    max_it( Parameter::getValue< int >( "fem.ode.maxiterations" , 16 ) ),
    sigma( Parameter::getValue< double >( "fem.ode.cflincrease" , 1.1 ) )
  {
  }
  virtual PARDG::IterativeLinearSolver *linearSolver(PARDG::Communicator & comm) const
  {
    int cycles = Parameter::getValue< int >( "fem.ode.gmrescycles" , 15 );
    PARDG::IterativeLinearSolver* solver = new PARDG::GMRES(comm,cycles);
    double tol = Parameter::getValue< double >( "fem.ode.solver.tolerance" , 1e-6 );
    static const std::string errorTypeTable[]
      = { "absolute", "relative" };
    int errorType = Parameter::getEnum( "fem.ode.solver.errormeassure", errorTypeTable, 0 );
    solver->set_tolerance(tol,(errorType==1));
    int maxIter = Parameter::getValue< int >( "fem.ode.solver.iterations" , 1000 );
    solver->set_max_number_of_iterations(maxIter);
    return solver;
  }
  virtual double tolerance() const
  {
    return Parameter::getValue< double >( "fem.ode.tolerance" , 1e-8 );
  }
  virtual int iterations() const
  {
    return Parameter::getValue< int >( "fem.ode.iterations" , 1000 );
  }
  virtual int verbose() const
  {
    static const std::string verboseTypeTable[]
      = { "none", "cfl", "full" };
    return Parameter::getEnum( "fem.ode.verbose" , verboseTypeTable, 0 );
  }
  virtual bool cflFactor( const PARDG::ODESolver &ode,
                          const PARDG::IterativeLinearSolver &solver,
                          bool converged,
                          double &factor) const
  {
    const int iter = solver.number_of_iterations();
    factor = 1.;
    bool changed = false;
    if (converged) 
    {
      if (iter < min_it) 
      {
        factor = sigma;
        changed = true;
      }
      else if (iter > max_it) 
      {
        factor = (double)max_it/(sigma*(double)iter);
        changed = true;
      }
    }
    else
    {
      factor = 0.5;
      changed = true;
    }
    return changed;
  }
  const int min_it,max_it;
  const double sigma;
};
#ifdef USE_PARDG_ODE_SOLVER
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

  //! return reference to real operator 
  const Operator& op() const { return op_; }

protected:
  // operator to call 
  Operator& op_;
  // discrete function space 
  const SpaceType& space_;
};


/**
   @ingroup ODESolver
   @{
 **/

/* \brief Base class for ParDG ode solvers */
template<class DestinationImp>
class ParDGOdeSolverBase : public OdeSolverInterface<DestinationImp>
{
  ParDGOdeSolverBase(const ParDGOdeSolverBase& );
public:
  //! type of destination function 
  typedef DestinationImp DestinationType; 

  //! type of discretization operator 
  typedef SpaceOperatorInterface<DestinationType> OperatorType;

protected:
  //! constructor
  ParDGOdeSolverBase(Dune::TimeProviderBase& tp, 
                     const int order ) :
    timeProvider_( tp ),
    comm_(PARDG::Communicator::instance()),
    order_( order ),
    initialized_( false ),
    odeSolver_( 0 )
  {
  }
 
public:  
  // initialize time step size 
  void initialize ( const DestinationType& U0 ) 
  {
    // initialized dt on first call
    if ( ! initialized_ )     
    {
      // create instance of ode solver 
      if( ! odeSolver_ ) odeSolver_ = createOdeSolver();

      // create temporary function
      DestinationType tmp( U0 );
      // apply operator once 
      spaceOperator()( U0, tmp );
      // set initial time step 
      timeProvider_.provideTimeStepEstimate( spaceOperator().timeStepEstimate() );
      // now it's initialized 
      initialized_ = true;
    }
  }

protected:  
  //! destructor  
  ~ParDGOdeSolverBase() { delete odeSolver_; odeSolver_ = 0; }

  //! return reference to ode solver 
  PARDG::ODESolver& odeSolver() 
  {
    assert( odeSolver_ );
    return *odeSolver_;
  }
  
  void printmyInfo(string filename) const 
  {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "Explicit ODE solver, steps: " << order_ << "\n\n";
    ofs.close();
    //this->op_.printmyInfo(filename);
  }

protected:  
  //! return ode solver object 
  virtual PARDG::ODESolver* createOdeSolver() = 0;

  //! return reference to discretization operator for time step initialization 
  virtual const OperatorType& spaceOperator() const = 0;
  
protected:
  Dune::TimeProviderBase& timeProvider_;
  PARDG::Communicator & comm_;
  const int order_;
  bool initialized_;

private:
  // use method odeSolver for access 
  PARDG::ODESolver* odeSolver_;
};

#if 0
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
    ode_(0),
    initialized_(false)
  {
    switch (pord) {
      case 1: ode_ = new PARDG::ExplicitEuler(comm_,expl_); break;
      case 2: ode_ = new PARDG::ExplicitTVD2(comm_,expl_); break;
      case 3: ode_ = new PARDG::ExplicitTVD3(comm_,expl_); break;
      case 4: ode_ = new PARDG::ExplicitRK4(comm_,expl_); break;
      default : ode_ = new PARDG::ExplicitBulirschStoer(comm_,expl_,7);
                std::cerr << "Runge-Kutta method of this order not implemented.\n" 
                          << "Using 7-stage Bulirsch-Stoer scheme.\n"
                          << std::endl;
    }

    if(verbose)
    {
      ode_->DynamicalObject::set_output(cout);
    }
  }
 
  // initialize time step size 
  bool initialize ( const DestinationType& U0 ) 
  {
    // initialized dt on first call
    if ( ! initialized_ )     
    {
      DestinationType tmp(U0);
      op_( U0, tmp );
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
  PARDG::ODESolver* ode_;
  bool initialized_;
};
#endif

//! ExplicitOdeSolver 
template<class DestinationImp>
class ExplicitOdeSolver : 
  public ParDGOdeSolverBase< DestinationImp >
{
  typedef ParDGOdeSolverBase< DestinationImp > BaseType; 

protected:  
  using BaseType :: order_;
  using BaseType :: comm_;
  using BaseType :: initialized_;
  using BaseType :: timeProvider_;
  using BaseType :: odeSolver ;
  using BaseType :: initialize ;

public:
  typedef typename BaseType :: OperatorType    OperatorType;
  typedef typename BaseType :: DestinationType DestinationType; 

  //! constructor 
  ExplicitOdeSolver(OperatorType& op, Dune :: TimeProviderBase &tp, const int order, bool verbose = false) :
    BaseType(tp, order ),
    expl_( op ),
    verbose_( verbose )
  {}

  //! destructor 
  virtual ~ExplicitOdeSolver() {}
 
  //! solve system 
  void solve( DestinationType& U0 ) 
  {
    // initialize 
    if( ! initialized_ ) 
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
    const bool convergence = odeSolver().step(time, dt , u);

    // set time step estimate of operator 
    timeProvider_.provideTimeStepEstimate( expl_.op().timeStepEstimate() );
    
    assert( convergence );
    if(! convergence ) 
    {
      timeProvider_.invalidateTimeStep();
      std::cerr << "No Convergence of ExplicitOdeSolver! \n";
    }
  }

protected:
  //! for initialization 
  const OperatorType& spaceOperator() const { return expl_.op(); }

  //! create explicit ode solver 
  PARDG::ODESolver* createOdeSolver() 
  {
    PARDG::ODESolver* odeSolver = 0;
    switch ( order_ ) 
    {
      case 1  : odeSolver = new PARDG::ExplicitEuler(comm_, expl_); break;
      case 2  : odeSolver = new PARDG::ExplicitTVD2 (comm_, expl_); break;
      case 3  : odeSolver = new PARDG::ExplicitTVD3 (comm_, expl_); break;
      case 4  : odeSolver = new PARDG::ExplicitRK4  (comm_, expl_); break;
      default : odeSolver = new PARDG::ExplicitBulirschStoer(comm_, expl_, 7);
                std::cerr << "Runge-Kutta method of this order not implemented.\n" 
                          << "Using 7-stage Bulirsch-Stoer scheme.\n"
                          << std::endl;
    }

    if( verbose_ )
    {
      odeSolver->DynamicalObject::set_output(cout);
    }
    assert( odeSolver );
    return odeSolver;
  }
 
protected:
  OperatorWrapper<OperatorType> expl_;
  const bool verbose_;
};

///////////////////////////////////////////////////////
//
//  --ImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class ImplicitOdeSolver : 
  public ParDGOdeSolverBase< DestinationImp >
{
  typedef ParDGOdeSolverBase< DestinationImp > BaseType; 
protected:
  using BaseType :: order_;
  using BaseType :: comm_;
  using BaseType :: initialized_;
  using BaseType :: timeProvider_;
  using BaseType :: odeSolver ;
  using BaseType :: initialize ;

public:
  typedef typename BaseType :: OperatorType    OperatorType;
  typedef typename BaseType :: DestinationType DestinationType; 

  ImplicitOdeSolver(OperatorType& op, 
                    Dune::TimeProviderBase& tp,
                    const int order, 
                    bool verbose ) DUNE_DEPRECATED :
    BaseType(tp, order),
    impl_( op ),
    linsolver_( 0 ),
    param_( ODEParameters().clone() ),
    verbose_( param_->verbose() ),
    cfl_(1.0)
  {
  }

  ImplicitOdeSolver(OperatorType& op, 
                    Dune::TimeProviderBase& tp,
                    const int order,
                    const ODEParameters& parameter= ODEParameters() ) :
    BaseType(tp, order),
    impl_( op ),
    linsolver_( 0 ),
    param_( parameter.clone() ),
    verbose_( parameter.verbose() ),
    cfl_(1.0)
  {
  }

protected:  
  //! return reference to parameter class 
  const ODEParameters& parameter() const { assert( param_ ); return *param_; }

  //! for initialization 
  const OperatorType& spaceOperator() const { return impl_.op(); }

  //! create implicit ode solver     
  PARDG::ODESolver* createOdeSolver()
  {
    linsolver_ = parameter().linearSolver( comm_ );
    assert( linsolver_ );

    PARDG::DIRK* odeSolver = 0;
    switch ( order_ ) 
    {
      case 1: odeSolver = new PARDG::ImplicitEuler(comm_, impl_); break;
      case 2: odeSolver = new PARDG::Gauss2(comm_, impl_); break;
      case 3: odeSolver = new PARDG::DIRK3 (comm_, impl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }

    assert( odeSolver );

    odeSolver->set_linear_solver( *linsolver_ );
    odeSolver->set_tolerance( parameter().tolerance() );
    odeSolver->set_max_number_of_iterations( parameter().iterations() );
    
    if( verbose_ == 2 ) 
    {
      odeSolver->IterativeSolver::set_output(cout);
      odeSolver->DynamicalObject::set_output(cout);
    }
    return odeSolver;
  }
  
  virtual ~ImplicitOdeSolver() 
  {
    delete linsolver_; linsolver_ = 0;
    delete param_;     param_ = 0;
  }
  
  //! return reference to ode solver 
  PARDG::DIRK& implicitSolver() 
  {
    return static_cast<PARDG::DIRK&> (odeSolver());
  }
  
public:  
  //! solve 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ImplicitOdeSolver wasn't initialized before first call!");
    }

    const double dt   = timeProvider_.deltaT();
    assert( dt > 0.0 );
    const double time = timeProvider_.time();

     // get pointer to solution
    double* u = U0.leakPointer();
      
    const bool convergence = odeSolver().step(time , dt , u);

    double factor;
    bool changed = parameter().cflFactor(odeSolver(),
                                         *(linsolver_),
                                         convergence, factor);
    cfl_ *= factor;
    if (convergence)
    {
      timeProvider_.provideTimeStepEstimate( cfl_ * impl_.op().timeStepEstimate() );

      if( changed && verbose_ >= 1 )
        derr << " New cfl number is: "<< cfl_ << " (number of iterations ("
             << "linear: " << linsolver_->number_of_iterations() 
             << ", ode: " << implicitSolver().number_of_iterations()
             << ")"
             << std::endl;
    } 
    else 
    {
      timeProvider_.provideTimeStepEstimate( cfl_ * dt );
      timeProvider_.invalidateTimeStep();
     
      // output only in verbose mode 
      if( verbose_ >= 1 )
      {
        derr << "No convergence: New cfl number is "<< cfl_ << std :: endl;
      }
    }

    linsolver_->reset_number_of_iterations();
  }

protected:  
  OperatorWrapper<OperatorType> impl_;
  PARDG::IterativeLinearSolver* linsolver_;
  const ODEParameters* param_;
  const int verbose_;
  double cfl_;
}; // end ImplicitOdeSolver

//////////////////////////////////////////////////////////////////
//
//  --SemiImplTimeStepperBase
//
//////////////////////////////////////////////////////////////////
template<class OperatorExpl, class OperatorImpl>
class SemiImplTimeStepperBase
{
  typedef typename OperatorExpl :: DestinationType DestinationType; 
public:
  SemiImplTimeStepperBase(OperatorExpl& explOp, OperatorImpl & implOp, Dune :: TimeProviderBase &tp, 
                      int pord, 
                      const ODEParameters& parameter=ODEParameters()) :
    ord_(pord),
    comm_(PARDG::Communicator::instance()),
    explOp_(explOp),
    implOp_(implOp),
    expl_(explOp),
    impl_(implOp),
    ode_(0),
    linsolver_(0),
    initialized_(false),
    verbose_(parameter.verbose()),
    param_(parameter.clone())
  {
    linsolver_ = parameter.linearSolver( comm_ );
    switch (pord) {
      case 1: ode_=new PARDG::SemiImplicitEuler(comm_,impl_,expl_); break;
      case 2: ode_=new PARDG::IMEX_SSP222(comm_,impl_,expl_); break;
      case 3: ode_=new PARDG::SIRK33(comm_,impl_,expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }
    ode_->set_linear_solver(*linsolver_);
    ode_->set_tolerance( parameter.tolerance() );
    ode_->set_max_number_of_iterations( parameter.iterations() );
    
    if( verbose_ == 2 ) 
    {
      ode_->IterativeSolver::set_output(cout);
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
      this->explOp_(U0,tmp);
      initialized_ = true;
      return true;
    }
    return false;
  }
  //! destructor 
  ~SemiImplTimeStepperBase() {delete ode_;delete linsolver_;}
  
  // return reference to ode solver 
  PARDG::SIRK& odeSolver() 
  {
    assert( ode_ );
    return *ode_;
  }

  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "Semi Implicit ODE solver, steps: " << this->ord_ << "\n\n";
    ofs.close();
    // this->op_.printmyInfo(filename);
  }
  
protected:
  int ord_;
  PARDG::Communicator & comm_;   
  const OperatorExpl& explOp_;
  const OperatorImpl& implOp_;
  OperatorWrapper<OperatorExpl> expl_;
  OperatorWrapper<OperatorImpl> impl_;
  PARDG::SIRK* ode_;
  PARDG::IterativeLinearSolver* linsolver_;
  bool initialized_;
  int verbose_;
  const ODEParameters* param_;
};


///////////////////////////////////////////////////////
//
//  --SemiImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class SemiImplicitOdeSolver : 
  public OdeSolverInterface<DestinationImp> ,
  public SemiImplTimeStepperBase<SpaceOperatorInterface<DestinationImp>, SpaceOperatorInterface<DestinationImp> > 
{
  typedef DestinationImp DestinationType;
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef SemiImplTimeStepperBase<OperatorType, OperatorType> BaseType;

protected:
  using BaseType :: implOp_;
  using BaseType :: explOp_;
  using BaseType :: odeSolver ;
  using BaseType :: linsolver_ ;
  using BaseType :: verbose_ ;
  using BaseType :: initialized_ ;

  Dune :: TimeProviderBase &timeProvider_;
  double cfl_;

public:
  SemiImplicitOdeSolver(OperatorType& explOp, OperatorType& implOp, Dune::TimeProviderBase& tp,
                    int pord, bool verbose ) DUNE_DEPRECATED :
    BaseType(explOp,implOp,tp,pord),
    timeProvider_(tp),
    cfl_(1.0)
  {
  }
  SemiImplicitOdeSolver(OperatorType& explOp, OperatorType& implOp, Dune::TimeProviderBase& tp,
                    int pord,
                    const ODEParameters& parameter=ODEParameters()) :
    BaseType(explOp,implOp,tp,pord,parameter),
    timeProvider_(tp),
    cfl_(1.0)
  {
  }

  virtual ~SemiImplicitOdeSolver() {}
  
  void provideTimeStepEstimate() const 
  {
    const double dtEst = std::min( explOp_.timeStepEstimate(),
                                   implOp_.timeStepEstimate() ); 

    // set time step estimate of operator 
    timeProvider_.provideTimeStepEstimate( cfl_ * dtEst );
  }

  //! initialize solver 
  void initialize(const DestinationType& U0)
  {
    // initialize solver 
    BaseType :: initialize (U0);

    // set time step size 
    provideTimeStepEstimate();
  }

  //! solve 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ImplicitOdeSolver wasn't initialized before first call!");
    }
    
    const double dt   = timeProvider_.deltaT();
    assert( dt > 0.0 );
    const double time = timeProvider_.time();

     // get pointer to solution
    double* u = U0.leakPointer();
      
    const bool convergence = odeSolver().step(time , dt , u);

    double factor;
    bool changed = BaseType::param_->cflFactor(odeSolver(),
                                               *(linsolver_),
                                               convergence,factor);
    cfl_ *= factor;
    if (convergence)
    {
      // set time step size 
      provideTimeStepEstimate();

      if( changed && verbose_ >= 1 )
        derr << " New cfl number is: "<< cfl_ << " (number of iterations: "
             << linsolver_->number_of_iterations() 
             << std::endl;
    } 
    else 
    {
      // reset time step size 
      timeProvider_.provideTimeStepEstimate( cfl_ * dt );
      timeProvider_.invalidateTimeStep();
     
      // output only in verbose mode 
      if( verbose_ >= 1 )
      {
        derr << "No convergence: New cfl number is "<< cfl_ << std :: endl;
      }
    }
    linsolver_->reset_number_of_iterations();
  }

}; // end SemiImplicitOdeSolver

#endif // USE_PARDG_ODE_SOLVER
 
/**
 @} 
**/

} // end namespace DuneODE

#undef USE_EXTERNAL_BLAS
#endif
