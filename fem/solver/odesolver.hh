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
#include <dune/fem/misc/timeutility.hh>

//#if HAVE_BLAS 
//#define USE_EXTERNAL_BLAS
//#endif

#include "ode/blas.hpp"

namespace DuneODE {

#include "ode/communicator.hpp"    
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"

// use Dennis namespace pardg
using namespace pardg;

} // end namespace DuneODE

namespace DuneODE {
  using namespace Dune;
  using namespace std;

template <class Operator>
class OperatorWrapper : public Function 
{
  typedef typename Operator::DestinationType DestinationType;
 public:
  OperatorWrapper(const Operator& op, TimeProvider& tp) 
    : op_(op) , tp_(tp) {}

  //! apply operator 
  void operator()(const double *u, double *f, int i = 0) {
    // create fake argument 
    DestinationType arg("ARG",op_.space(),u);
    // create fake destination 
    DestinationType dest("DEST",op_.space(),f);
    
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
  // time provider 
  TimeProvider& tp_;
};
  /** @defgroup ODESolver ODE Solver
   *  @ingroup OperatorCommon
   @{
   **/

template <class Operator>
class OdeSolverInterface {
protected:
  OdeSolverInterface () {}    

  typedef typename Operator :: DestinationType DestinationType;
public:
  //! destructor 
  virtual ~OdeSolverInterface () {}
  
  //! solve system 
  virtual void solve(DestinationType&) = 0;
};

template<class Operator>
class ExplTimeStepperBase 
{
  typedef typename Operator::DestinationType DestinationType; 
public:
  ExplTimeStepperBase(Operator& op, TimeProvider& tp, 
                      int pord, bool verbose) :
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

  ~ExplTimeStepperBase() { delete ode_; }

  // return reference to ode solver 
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
template<class Operator>
class ExplicitOdeSolver : 
  public OdeSolverInterface<Operator> ,
  public ExplTimeStepperBase<Operator>  
{
  typedef ExplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType; 
 public:
  //! constructor 
  ExplicitOdeSolver(Operator& op, TimeProvider& tp, int pord, bool verbose = false) :
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
  
  //! solve system 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( BaseType::initialize(U0) ) return ;

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
template<class Operator>
class ImplicitOdeSolver : 
  public OdeSolverInterface<Operator> ,
  public ImplTimeStepperBase<Operator> 
{
  typedef ImplTimeStepperBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType;
public:
  ImplicitOdeSolver(Operator& op, TimeProvider& tp,
                    int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {
  }

  virtual ~ImplicitOdeSolver() {}
  
  void solve(DestinationType& U0) 
  {
    // for first call only calculate time step estimate 
    if( BaseType::initialize(U0) ) return ;

    bool convergence = false;
    int cycle = 0;
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
        derr << "New cfl number is "<< timeProvider_.cfl() << "\n";
      }

      ++cycle;
      if( cycle > 5 ) 
      {
        derr << "No Convergence of implicit ODE solver! \n";
        abort();
      }
    }
  }

private:
  TimeProvider& timeProvider_;
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

template<class Operator>
class ExplRungeKuttaBase 
{
public:
  typedef typename Operator::DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
private:
  std::vector< std::vector<double> > a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<DestinationType*> Upd;
protected:  
  const int ord_;

public:
  ExplRungeKuttaBase(Operator& op, TimeProvider& tp, 
                     int pord, bool verbose = true ) :
    a(0),b(0),c(0), Upd(0),
    ord_(pord),
    op_(op),
    tp_(tp),
    initialized_(false)
  {
    assert(ord_>0);
    a.resize(ord_);
    for (int i=0; i<ord_; ++i)
    {
      a[i].resize(ord_);
    }
    b.resize(ord_); 
    c.resize(ord_); 
    
    switch (ord_) {
    case 4 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;    a[0][3]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;    a[1][3]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;    a[2][3]=0.;
      a[3][0]=1./6.;  a[3][1]=1./6.;  a[3][2]=2./3.; a[3][3]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;    b[3]=0.;
      c[0]=0.;        c[1]=1.0;       c[2]=0.5;      c[3]=1.0;
      break;
    case 3 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;
      c[0]=0.;        c[1]=1;         c[2]=0.5;
      break;
    case 2 :
      a[0][0]=0.;     a[0][1]=0.;
      a[1][0]=1.0;    a[1][1]=0.;
      b[0]=0.5;       b[1]=0.5;
      c[0]=0;         c[1]=1;
      break;
    case 1:
      a[0][0]=0.;
      b[0]=1.;
      c[0]=0.;
      break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
      << std::endl;
              abort();
    }

    // create update memory 
    for (int i=0; i<ord_; ++i)
    {
      Upd.push_back(new DestinationType("URK",op_.space()) );
    }
    Upd.push_back(new DestinationType("Ustep",op_.space()) );
  }

  ~ExplRungeKuttaBase()
  {
    for(size_t i=0; i<Upd.size(); ++i) 
      delete Upd[i];
  }

  // apply operator once to get dt estimate 
  bool initialize(const DestinationType& U0)
  {
    if( ! initialized_ ) 
    {
      // Compute Steps
      op_(U0, *(Upd[0]));
      initialized_ = true;
      return true;
    }
    return false;
  }
  
  void solve(DestinationType& U0) 
  {
    // time might change 
    tp_.unlock();
    
    // get cfl * timeStepEstimate 
    const double dt = tp_.deltaT();
    // get time 
    const double t = tp_.time();

    // Compute Steps
    op_(U0, *(Upd[0]));
    
    for (int i=1; i<ord_; ++i) 
    {
      (Upd[ord_])->assign(U0);
      for (int j=0; j<i ; ++j) 
      {
        (Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      }

      // set new time 
      tp_.setTime( t + c[i]*dt );

      // apply operator 
      op_(*(Upd[ord_]),*(Upd[i]));
    }

    // Perform Update
    for (int j=0; j<ord_; ++j) 
    {
      U0.addScaled(*(Upd[j]),(b[j]*dt));
    }
    
    // restore global time 
    tp_.lock();
  }

protected:
  // operator to solve for 
  const Operator& op_;
  // time provider 
  TimeProvider& tp_;
  // init flag 
  bool initialized_;
};

template<class Operator>
class ExplRungeKutta : public TimeProvider , 
                       public ExplRungeKuttaBase<Operator> 
{
  typedef ExplRungeKuttaBase<Operator> BaseType;
public:
  typedef typename Operator :: DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 

public:
  ExplRungeKutta(Operator& op,int pord,double cfl, bool verbose = true ) :
    TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(op.space().grid().comm(),*this), 
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
  }

  void initialize(const DestinationType& U0)
  {
    if(! this->initialized_)
    {
      // initialize 
      BaseType :: initialize(U0);
    
      // global min of dt and reset of dtEstimate 
      this->tp_.syncTimeStep();
    }
  }
    
  double solve(typename Operator::DestinationType& U0) 
  {
    initialize( U0 );
    
    // solve ode 
    BaseType :: solve (U0);
    
    // calls setTime ( t + dt ); 
    this->tp_.augmentTime();
    
    // global min of dt and reset of dtEstimate 
    this->tp_.syncTimeStep();
    
    return this->tp_.time();
  }
  
  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,this->op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ExplRungeKutta, steps: " << this->ord_ << "\n\n";
    ofs << "                cfl: " << this->tp_.cfl() << "\\\\\n\n";
    ofs.close();
    this->op_.printmyInfo(filename);
  }

private:
  // TimeProvider with communicator 
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  double savetime_;
  int savestep_;
};

//! ExplicitRungeKuttaSolver 
template<class Operator>
class ExplicitRungeKuttaSolver : 
  public OdeSolverInterface<Operator> ,
  public ExplRungeKuttaBase<Operator>  
{
  typedef ExplRungeKuttaBase<Operator> BaseType;
  typedef typename Operator::DestinationType DestinationType; 
 public:
  //! constructor 
  ExplicitRungeKuttaSolver(Operator& op, TimeProvider& tp, int pord, bool verbose = false) :
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
  virtual ~ExplicitRungeKuttaSolver() {}
  
  //! solve system 
  void solve(DestinationType& U0) 
  {
    // on first call just apply operator once 
    if( BaseType :: initialize(U0) ) return ;

    // solve ode 
    BaseType :: solve(U0);
  }

private:
  TimeProvider& timeProvider_;
};


//////////////////////////////////////////////////////////
//
// Operator Interface to use linear solvers from DuneODE
//
//////////////////////////////////////////////////////////
template <class OperatorImp>
class SolverInterfaceImpl : public Function 
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
/** @} **/
}

#undef USE_EXTERNAL_BLAS
#endif
