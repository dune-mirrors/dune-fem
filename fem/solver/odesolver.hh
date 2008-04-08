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
#include <dune/fem/misc/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

//- include runge kutta ode solver 
#include <dune/fem/solver/rungekutta.hh>

// include headers of pardg 
#include "pardg.hh"

namespace DuneODE {

#ifdef USE_PARDG_ODE_SOLVER

template <class Operator>
class OperatorWrapper : public pardg::Function 
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
                      Dune::TimeProvider& tp, 
                      int pord, 
                      bool verbose) :
    ord_(pord),
    comm_(pardg::Communicator::instance()),
    op_(op),
    expl_(op),
    ode_(0),
    initialized_(false)
  {
    switch (pord) {
      case 1: ode_ = new pardg::ExplicitEuler(comm_,expl_); break;
      case 2: ode_ = new pardg::ExplicitTVD2(comm_,expl_); break;
      case 3: ode_ = new pardg::ExplicitTVD3(comm_,expl_); break;
      case 4: ode_ = new pardg::ExplicitRK4(comm_,expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }

    if(verbose)
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
  pardg::ODESolver& odeSolver() 
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
  pardg::Communicator & comm_;
  const Operator& op_;
  OperatorWrapper<Operator> expl_;
  pardg::ODESolver* ode_;
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
  ExplicitOdeSolver(OperatorType& op, Dune::TimeProvider& tp, int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {
    if(verbose) 
    {
      std::cout << "ExplicitOdeSolver: cfl = " << tp.cfl() << "!\n";
    } 
  }

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

private:
  Dune::TimeProvider& timeProvider_;
};

template<class Operator>
class ExplTimeStepper : public Dune::TimeProvider, 
                        public ExplTimeStepperBase<Operator>  
{
  typedef ExplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType; 
  typedef typename DestinationType :: DiscreteFunctionSpaceType 
    :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
public:
  ExplTimeStepper(Operator& op,int pord, double cfl) :
    Dune::TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,false),
    tp_(this->op_.space().grid().comm(), *this ),
    savestep_(1),
    savetime_(0.0)
  {
    assert( this->cfl() <= 1.0 );
  }
  
  ExplTimeStepper(Operator& op,int pord, double cfl, bool verbose ) :
    Dune::TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,true),
    tp_(this->op_.space().grid().comm(), *this ),
    savestep_(1),
    savetime_(0.0)
  {
    assert( this->cfl() <= 1.0 );
  }

  ExplTimeStepper(Operator& op,int pord, double cfl, double startTime) :
    Dune::TimeProvider(startTime,cfl),
    BaseType(op,*this,pord,false),
    tp_(this->op_.space().grid().comm(), *this ),
    savestep_(1),
    savetime_(startTime)
  {
    assert( this->cfl() <= 1.0 );
  }

  // initialize 
  void initialize(const DestinationType& U0)
  {
    if( ! this->initialized_ ) 
    {
      BaseType::initialize(U0);

      // set time step estimate of operator 
      tp_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
    
      // global min of dt 
      tp_.syncTimeStep(); 

      this->initialized_ = true;
    }
  }
 
  // solve ode 
  double solve(typename Operator::DestinationType& U0) 
  {
    // initialize time step size 
    initialize(U0);
    
    // reset estimate 
    tp_.resetTimeStepEstimate();
    
    // get dof array 
    double* u=U0.leakPointer();

    assert( tp_.deltaT() > 0.0 );

    const double t  = tp_.time();
    const double dt = tp_.deltaT();

    // solve ode 
    const bool convergence = this->odeSolver().step(t, dt, u);

    // set time step estimate of operator 
    tp_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
    
    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ExplTimeStepper! \n";
      tp_.invalidateTimeStep();
    }

    // increase time step 
    tp_.next();

    // return time 
    return tp_.time();
  }

  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) 
    {
      printSGrid(time(),savestep_*10+nr,this->op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
 private:
  Dune::ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
};

//////////////////////////////////////////////////////////////////
//
//  --ImplTimeStepperBase
//
//////////////////////////////////////////////////////////////////
template<class Operator>
class ImplTimeStepperBase
{
  typedef typename Operator :: DestinationType DestinationType; 
public:
  ImplTimeStepperBase(Operator& op, Dune::TimeProvider& tp, 
                      int pord, bool verbose) :
    ord_(pord),
    comm_(pardg::Communicator::instance()),
    op_(op),
    impl_(op),
    ode_(0),
    linsolver_(comm_,cycle),
    initialized_(false)
  {
    linsolver_.set_tolerance(1.0e-8,false);
    linsolver_.set_max_number_of_iterations(10000);
    switch (pord) 
    {
      case 1: ode_ = new pardg::ImplicitEuler(comm_,impl_); break;
      case 2: ode_ = new pardg::Gauss2(comm_,impl_); break;
      case 3: ode_ = new pardg::DIRK3(comm_,impl_); break;
      //case 4: ode_ = new pardg::ExplicitRK4(comm,expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }
    ode_->set_linear_solver(linsolver_);
    ode_->set_tolerance(1.0e-6);
    ode_->set_max_number_of_iterations(15);
    if( verbose ) 
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
      this->op_(U0,tmp);
      initialized_ = true;
      return true;
    }
    return false;
  }
  //! destructor 
  ~ImplTimeStepperBase() {delete ode_;}
  
  // return reference to ode solver 
  pardg::DIRK& odeSolver() 
  {
    assert( ode_ );
    return *ode_;
  }
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "Implicit ODE solver, steps: " << this->ord_ << "\n\n";
    ofs.close();
    // this->op_.printmyInfo(filename);
  }
  
protected:
  int ord_;
  pardg::Communicator & comm_;   
  const Operator& op_;
  OperatorWrapper<Operator> impl_;
  pardg::DIRK* ode_;
  pardg::GMRES linsolver_;
  enum { cycle = 20 };
  bool initialized_;
};


///////////////////////////////////////////////////////
//
//  --ImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class ImplicitOdeSolver : 
  public OdeSolverInterface<DestinationImp> ,
  public ImplTimeStepperBase<SpaceOperatorInterface<DestinationImp> > 
{
  typedef DestinationImp DestinationType;
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef ImplTimeStepperBase<OperatorType> BaseType;
public:
  ImplicitOdeSolver(OperatorType& op, Dune::TimeProvider& tp,
                    int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp),
    cfl_(1.0)
  {
  }

  virtual ~ImplicitOdeSolver() {}
  
  //! initialize solver 
  void initialize(const DestinationType& U0)
  {
    // get current cfl estimate 
    cfl_ = timeProvider_.cfl();
    // initialize solver 
    BaseType :: initialize (U0);

    // set time step estimate of operator 
    timeProvider_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
  }

  //! solve 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! this->initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ImplicitOdeSolver wasn't initialized before first call!");
    }

    bool convergence = false;
    int cycle = 0;

    // if cfl has been changed 
    // then try to reach old value 
    if( cfl_ > timeProvider_.cfl() ) 
    {
      double cfl = 2.0 * timeProvider_.cfl();
      timeProvider_.setCfl(cfl);
    }
    
    while( !convergence )
    {
      const double dt   = timeProvider_.deltaT();
      assert( dt > 0.0 );
      const double time = timeProvider_.time();

      // get pointer to solution
      double* u = U0.leakPointer();
      
      convergence = this->odeSolver().step(time , dt , u);
      // set time step estimate of operator 
      timeProvider_.provideTimeStepEstimate( this->op_.timeStepEstimate() );
      
      if(!convergence) 
      {
        double cfl = 0.5 * timeProvider_.cfl();
        timeProvider_.setCfl(cfl); 
        // output only on rank 0
        if(U0.space().grid().comm().rank() == 0 )
        {
          derr << "New cfl number is "<< timeProvider_.cfl() << "\n";
        }
      }

      ++cycle;
      if( cycle > 25 ) 
      {
        DUNE_THROW(InvalidStateException,"ImplicitOdeSolver: no convergence of solver!");
      }
    }
    if (this->odeSolver().number_of_iterations()<10) {
      double cfl = 2.0 * timeProvider_.cfl();
      timeProvider_.setCfl(cfl); 
      if(U0.space().grid().comm().rank() == 0 )
      {
        derr << "New cfl number is "<< timeProvider_.cfl() << "\n";
      }
    }

  }

private:
  Dune::TimeProvider& timeProvider_;
  double cfl_;
};


///////////////////////////////////////////////////////
//
//  --ImplTimeStepper 
//
///////////////////////////////////////////////////////
template<class Operator>
class ImplTimeStepper : public Dune::TimeProvider ,
                        public ImplTimeStepperBase<Operator> 
{
  typedef ImplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType;
  typedef typename  DestinationType ::
      DiscreteFunctionSpaceType :: GridType :: Traits ::
      CollectiveCommunication DuneCommunicatorType; 
public:
  ImplTimeStepper(Operator& op,int pord,double cfl, bool verbose = false) :
    Dune::TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(this->op_.space().grid().comm(),*this),
    savestep_(1),
    savetime_(0.0)
  {
  }
  
  // initialize 
  void initialize(const DestinationType& U0)
  {
    if( ! this->initialized_ ) 
    {
      // call initialize 
      BaseType::initialize(U0);

      // global min of dt 
      tp_.syncTimeStep(); 
      this->initialized_ = true;
    }
  }
  
  // solve ode 
  double solve(DestinationType& U0) 
  {
    // initialize 
    initialize(U0);
    
    // reset dt estimate
    tp_.resetTimeStepEstimate();

    // get pointer to dof array 
    double* u=U0.leakPointer();
    
    assert( tp_.deltaT() > 0.0 );
    
    const double t  = tp_.time();
    const double dt = tp_.deltaT();

    // call ode solver  
    const bool convergence =  this->odeSolver().step(t, dt, u);

    // set time step estimate of operator 
    tp_.provideTimeStepEstimate( this->op_.timeStepEstimate() );

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ImplTimeStepper! \n";
      abort();
    }

    // next time step 
    tp_.next();

    return tp_.time();
  }

  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,this->op_.space(), U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  
private:
  Dune::ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
};


//////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////
template<class OperatorExpl,class OperatorImpl>
class SemiImplTimeStepper : public Dune::TimeProvider 
{
  typedef OperatorExpl Operator;
  typedef typename  Operator :: DestinationType ::
      DiscreteFunctionSpaceType :: GridType :: Traits ::
      CollectiveCommunication DuneCommunicatorType; 
 public:
  SemiImplTimeStepper(OperatorExpl& op_expl,OperatorImpl& op_impl,
          int pord,double cfl, bool verbose = false ) :
    Dune::TimeProvider(0.0,cfl),
    ord_(pord),
    comm_(pardg::Communicator::instance()),
    opexpl_(op_expl),
    opimpl_(op_impl),
    expl_(op_expl),
    impl_(op_impl),
    ode_(0),
    // linsolver_(comm_,cycle),
    linsolver_(comm_),
    tp_(this->opexpl_.space().grid().comm(),*this),
    savestep_(1),
    savetime_(0.0), 
    initialized_(false)
  {
    // linsolver_.set_tolerance(1.0e-8,false);
    linsolver_.set_tolerance(1.0e-3,false);
    linsolver_.set_max_number_of_iterations(1000);
    switch (pord) {
      case 1: ode_=new pardg::SemiImplicitEuler(comm_,impl_,expl_); break;
      case 2: ode_=new pardg::IMEX_SSP222(comm_,impl_,expl_); break;
      case 3: ode_=new pardg::SIRK33(comm_,impl_,expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }
    ode_->set_linear_solver(linsolver_);
    // ode_->set_tolerance(1.0e-6);
    ode_->set_tolerance(1.0e-2);
    if( verbose ) 
    { 
      ode_->IterativeSolver::set_output(cout);
      ode_->DynamicalObject::set_output(cout);
    }
  }
  
  ~SemiImplTimeStepper() {delete ode_;}

  double solve(typename Operator::DestinationType& U0) 
  {
    typedef typename Operator:: DestinationType :: DiscreteFunctionSpaceType :: 
          GridType :: Traits ::  CollectiveCommunication DuneCommunicatorType; 
    // const DuneCommunicatorType & duneComm = opexpl_.space().grid().comm();

    if ( ! initialized_ ) 
    {
      typename OperatorExpl::DestinationType tmp(U0);
      opexpl_(U0,tmp);
      
      // calculate global min of dt 
      tp_.syncTimeStep();

      initialized_ = true;
    }
    
    tp_.resetTimeStepEstimate();

    double* u=U0.leakPointer();
    
    const double t  = tp_.time();
    const double dt = tp_.deltaT();
    
    const bool convergence = ode_->step( t, dt, u);

    // set time step estimate of operator 
    tp_.provideTimeStepEstimate( this->opexpl_.timeStepEstimate() );

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ImplTimeStepper! \n";
      abort();
    }
    
    // next time step 
    tp_.next();
    
    return tp_.time();
  }

  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,opexpl_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  
  void printmyInfo(string filename) const 
  {
    std::ostringstream filestream;
    filestream << filename;
    {
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "SemiImplTimeStepper, steps: " << ord_ << "\n\n";
      ofs << "Explicite Operator:\\\\\n\n";
      ofs.close();
      // opexpl_.printmyInfo(filename);
    }
    {
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Implicite Operator:\\\\\n\n";
      ofs.close();
      // opimpl_.printmyInfo(filename);
    }
  }
 private:
  int ord_;
  pardg::Communicator & comm_;   
  const OperatorExpl& opexpl_;
  const OperatorImpl& opimpl_;
  OperatorWrapper<OperatorExpl> expl_;
  OperatorWrapper<OperatorImpl> impl_;
  pardg::SIRK* ode_;
  pardg::CG linsolver_;
  // pardg::GMRES linsolver_;
  enum { cycle = 20 };
  // Dune::TimeProvider with communicator 
  Dune::ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
  bool initialized_;
};
#endif
/**
 @} 
**/
}

#undef USE_EXTERNAL_BLAS
#endif
