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

#define USE_DENNIS_ODE_SOLVER

// if the preprocessor variable is defined, the ODE Solver from Dennis
// are used.
#ifdef USE_DENNIS_ODE_SOLVER

//#if HAVE_BLAS 
//#define USE_EXTERNAL_BLAS
//#endif

#include "ode/blas.hpp"

namespace DuneODE {

#include "ode/communicator.hpp"    
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"
#include "ode/bicgstab.hpp"

// use Dennis namespace pardg
using namespace pardg;

} // end namespace DuneODE
#endif

namespace DuneODE {
  using namespace Dune;
  using namespace std;

#ifdef USE_DENNIS_ODE_SOLVER

template <class Operator>
class OperatorWrapper : public Function 
{
  // type of discrete function 
  typedef typename Operator::DestinationType DestinationType;
  // type of discrete function space 
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
 public:
  //! constructor 
  OperatorWrapper(const Operator& op, TimeProvider& tp) 
    : op_(op) , space_(op_.space()) , tp_(tp) 
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
    tp_.setTime(time());
    
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
  const Operator& op_;
  // discrete function space 
  const SpaceType& space_;
  // time provider 
  TimeProvider& tp_;
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
                      TimeProvider& tp, 
                      int pord, 
                      bool verbose) :
    ord_(pord),
    comm_(Communicator::instance()),
    op_(op),
    expl_(op,tp),
    ode_(0),
    initialized_(false)
  {
    switch (pord) {
    case 1: ode_=new ExplicitEuler(comm_,expl_); break;
    case 2: ode_=new ExplicitTVD2(comm_,expl_); break;
    case 3: ode_=new ExplicitTVD3(comm_,expl_); break;
    case 4: ode_=new ExplicitRK4(comm_,expl_); break;
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
      DestinationType tmp("tmp",this->op_.space());
      this->op_(U0,tmp);
      initialized_ = true;
      return true;
    }
    return false;
  }

  //! destructor  
  ~ExplTimeStepperBase() { delete ode_; }

  //! return reference to ode solver 
  DuneODE::ODESolver& odeSolver() 
  {
    assert( ode_ );
    return *ode_;
  }
  
protected:
  int ord_;
  Communicator & comm_;
  const Operator& op_;
  OperatorWrapper<Operator> expl_;
  DuneODE::ODESolver* ode_;
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
  ExplicitOdeSolver(OperatorType& op, TimeProvider& tp, int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {
    // CFL upper estimate 
    double cfl = 0.45 / (2.0 * pord+1);

    // maximal allowed cfl number 
    tp.provideCflEstimate(cfl); 
    assert( tp.cfl() <= 1.0 );

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
    
    // time can is changeable now 
    timeProvider_.unlock();
    
    // call ode solver 
    const bool convergence = this->odeSolver().step(time, dt , u);

    // restore saved time 
    timeProvider_.lock();
    
    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ExplTimeStepper! \n";
      abort();
    }
  }

private:
  TimeProvider& timeProvider_;
};

template<class Operator>
class ExplTimeStepper : public TimeProvider, 
                        public ExplTimeStepperBase<Operator>  
{
  typedef ExplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType; 
  typedef typename DestinationType :: DiscreteFunctionSpaceType 
    :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
public:
  ExplTimeStepper(Operator& op,int pord, double cfl, bool verbose = false) :
    TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(this->op_.space().grid().comm(), *this ),
    savestep_(1),
    savetime_(0.0)
  {
    op.timeProvider(this);
    assert( this->cfl() <= 1.0 );
  }

  // initialize 
  void initialize(const DestinationType& U0)
  {
    if( ! this->initialized_ ) 
    {
      BaseType::initialize(U0);

      // global min of dt 
      tp_.syncTimeStep(); 
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

    // time can is changeable now 
    tp_.unlock();
    
    // solve ode 
    const bool convergence = this->odeSolver().step(t, dt, u);

    // restore saved time 
    tp_.lock();
    
    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ExplTimeStepper! \n";
      abort();
    }

    // set new time to t + deltaT 
    tp_.augmentTime();
    // check time increment 
    assert( std::abs( tp_.time() - (t+dt) ) <1e-10);
    
    // global min of new dt 
    tp_.syncTimeStep();
    
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
  void printmyInfo(string filename) const 
  {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ExplTimeStepper, steps: " << this->ord_ << "\n\n";
    ofs << "                 cfl: " << tp_.cfl() << "\\\\\n\n";
    ofs.close();
    this->op_.printmyInfo(filename);
  }
 private:
  ParallelTimeProvider<DuneCommunicatorType> tp_;
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
  ImplTimeStepperBase(Operator& op, TimeProvider& tp, 
                      int pord, bool verbose) :
    ord_(pord),
    comm_(Communicator::instance()),
    op_(op),
    impl_(op,tp),
    ode_(0),
    linsolver_(comm_,cycle),
    initialized_(false)
  {
    linsolver_.set_tolerance(1.0e-8,false);
    linsolver_.set_max_number_of_iterations(10000);
    switch (pord) {
    case 1: ode_=new ImplicitEuler(comm_,impl_); break;
    case 2: ode_=new Gauss2(comm_,impl_); break;
    case 3: ode_=new DIRK3(comm_,impl_); break;
      //case 4: ode_=new ExplicitRK4(comm,expl_); break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
      << std::endl;
      abort();
    }
    ode_->set_linear_solver(linsolver_);
    ode_->set_tolerance(1.0e-6);
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
      DestinationType tmp("tmp",this->op_.space());
      this->op_(U0,tmp);
      initialized_ = true;
      return true;
    }
    return false;
  }
  //! destructor 
  ~ImplTimeStepperBase() {delete ode_;}
  
  // return reference to ode solver 
  DuneODE::DIRK& odeSolver() 
  {
    assert( ode_ );
    return *ode_;
  }
  
protected:
  int ord_;
  Communicator & comm_;   
  const Operator& op_;
  OperatorWrapper<Operator> impl_;
  DuneODE::DIRK* ode_;
  DuneODE::GMRES linsolver_;
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
  ImplicitOdeSolver(OperatorType& op, TimeProvider& tp,
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
      
      // restore saved time 
      timeProvider_.unlock();
    
      convergence = this->odeSolver().step(time , dt , u);

      // restore saved time 
      timeProvider_.lock();
    
      if(!convergence) 
      {
        double cfl = 0.5 * timeProvider_.cfl();
        timeProvider_.provideCflEstimate(cfl); 
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
  }

private:
  TimeProvider& timeProvider_;
  double cfl_;
};


///////////////////////////////////////////////////////
//
//  --ImplTimeStepper 
//
///////////////////////////////////////////////////////
template<class Operator>
class ImplTimeStepper : public TimeProvider ,
                        public ImplTimeStepperBase<Operator> 
{
  typedef ImplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType;
  typedef typename  DestinationType ::
      DiscreteFunctionSpaceType :: GridType :: Traits ::
      CollectiveCommunication DuneCommunicatorType; 
public:
  ImplTimeStepper(Operator& op,int pord,double cfl, bool verbose = false) :
    TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(this->op_.space().grid().comm(),*this),
    savestep_(1),
    savetime_(0.0)
  {
    op.timeProvider(this);
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
    // time can is changeable now 
    tp_.unlock();

    // call ode solver  
    const bool convergence =  this->odeSolver().step(t, dt, u);

    // restore global time 
    tp_.lock();
    
    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ImplTimeStepper! \n";
      abort();
    }

    // calls setTime(t+ this->dt_);
    tp_.augmentTime(); 
    // check time increment 
    assert( std::abs( tp_.time() - (t+dt) ) <1e-10);
    
    // calculate global min of dt 
    tp_.syncTimeStep();

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
  
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ImplTimeStepper, steps: " << this->ord_ << "\n\n";
    ofs << "                 cfl: " << tp_.cfl() << "\\\\\n\n";
    ofs.close();
    this->op_.printmyInfo(filename);
  }
private:
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
};


//////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////
template<class OperatorExpl,class OperatorImpl>
class SemiImplTimeStepper : public TimeProvider 
{
  typedef OperatorExpl Operator;
  typedef typename  Operator :: DestinationType ::
      DiscreteFunctionSpaceType :: GridType :: Traits ::
      CollectiveCommunication DuneCommunicatorType; 
 public:
  SemiImplTimeStepper(OperatorExpl& op_expl,OperatorImpl& op_impl,
          int pord,double cfl, bool verbose = false ) :
    TimeProvider(0.0,cfl),
    ord_(pord),
    comm_(Communicator::instance()),
    opexpl_(op_expl),
    opimpl_(op_impl),
    expl_(op_expl,*this),
    impl_(op_impl,*this),
    ode_(0),
    linsolver_(comm_,cycle),
    // linsolver_(comm_),
    tp_(this->opexpl_.space().grid().comm(),*this),
    savestep_(1),
    savetime_(0.0), 
    initialized_(false)
  {
    op_expl.timeProvider(this);
    linsolver_.set_tolerance(1.0e-8,false);
    linsolver_.set_max_number_of_iterations(1000);
    switch (pord) {
    case 1: ode_=new SemiImplicitEuler(comm_,impl_,expl_); break;
    case 2: ode_=new IMEX_SSP222(comm_,impl_,expl_); break;
    case 3: ode_=new SIRK33(comm_,impl_,expl_); break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
      << std::endl;
      abort();
    }
    ode_->set_linear_solver(linsolver_);
    ode_->set_tolerance(1.0e-6);
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
    const DuneCommunicatorType & duneComm = opexpl_.space().grid().comm();

    if ( ! initialized_ ) 
    {
      typename OperatorExpl::DestinationType tmp("tmp",opexpl_.space());
      opexpl_(U0,tmp);
      
      // calculate global min of dt 
      tp_.syncTimeStep();

      initialized_ = true;
    }
    
    tp_.resetTimeStepEstimate();

    double* u=U0.leakPointer();
    
    const double t  = tp_.time();
    const double dt = tp_.deltaT();
    // time can is changeable now 
    tp_.unlock();
    
    const bool convergence = ode_->step( t, dt, u);
    assert(convergence);

    // restore global time 
    tp_.lock();
    
    // calls setTime(t+dt_);
    tp_.augmentTime();
    // check time increment 
    assert( std::abs( tp_.time() - (t+dt) ) <1e-10);
    
    // calculate global min of dt 
    tp_.syncTimeStep();

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
  
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    {
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "SemiImplTimeStepper, steps: " << ord_ << "\n\n";
      ofs << "                     cfl: " << tp_.cfl() << "\\\\\n\n";
      ofs << "Explicite Operator:\\\\\n\n";
      ofs.close();
      opexpl_.printmyInfo(filename);
    }
    {
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs << "Implicite Operator:\\\\\n\n";
      ofs.close();
      opimpl_.printmyInfo(filename);
    }
  }
 private:
  int ord_;
  Communicator & comm_;   
  const OperatorExpl& opexpl_;
  const OperatorImpl& opimpl_;
  OperatorWrapper<OperatorExpl> expl_;
  OperatorWrapper<OperatorImpl> impl_;
  DuneODE::SIRK* ode_;
  // DuneODE::CG linsolver_;
  DuneODE::GMRES linsolver_;
  enum { cycle = 20 };
  // TimeProvider with communicator 
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
  bool initialized_;
};
#endif

//////////////////////////////////////////////////////////
//
// Operator Interface to use linear solvers from DuneODE
//
//////////////////////////////////////////////////////////
template <class OperatorImp>
class SolverInterfaceImpl 
#ifdef USE_DENNIS_ODE_SOLVER
: public Function 
#endif
{
  const OperatorImp & op_;
  int size_; 
public:
  SolverInterfaceImpl(const OperatorImp & op, int size = 0) 
    : op_(op), size_(size) 
  {}

  void setSize( int size ) { size_ = size; }

  void operator () (const double *arg, double * dest, int i = 0 ) 
  {
    op_.multOEM(arg,dest);
  }
  
  void mult(const double *arg, double * dest) const
  {
    op_.multOEM(arg,dest);
  }
  
  int dim_of_argument(int i = 0) const 
  { 
    assert( i == 0 );
    return size_;
  }
  int dim_of_value(int i = 0) const 
  { 
    assert( i == 0 );
    return size_;
  }
};

/**
 @} 
**/
}

#undef USE_EXTERNAL_BLAS
#endif
