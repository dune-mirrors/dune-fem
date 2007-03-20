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

// include first, because of typedef 
#if HAVE_MPI 
#include "ode/communicator.cpp"    
#include "ode/buffer.cpp"
#else 
#include "ode/emptycommunicator.hpp"
#endif

  
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"
#include "ode/bulirsch_stoer.cpp"  
#include "ode/iterative_solver.cpp"  
#include "ode/ode_solver.cpp"     
#include "ode/sirk.cpp"   

#include "ode/matrix.cpp"            
#include "ode/qr_solver.cpp"   
#include "ode/ssp.cpp"
#include "ode/dirk.cpp"     
#include "ode/runge_kutta.cpp"
#include "ode/gmres.cpp"
#include "ode/fgmres.cpp"
#include "ode/bicgstab.cpp"
#include "ode/cg.cpp"
#include "ode/vector.cpp"

// use Dennis namespace pardg
using namespace pardg;
}

namespace DuneODE {
  using namespace Dune;
  using namespace std;

template <class Operator>
class OperatorWrapper : public Function {
 public:
  OperatorWrapper(const Operator& op) : op_(op) {}

  //! apply operator 
  void operator()(const double *u, double *f, int i = 0) {
    typename Operator::DestinationType arg("ARG",op_.space(),u);
    typename Operator::DestinationType dest("DEST",op_.space(),f);
    op_.setTime(time());
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
  const Operator& op_;
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
public:
  ExplTimeStepperBase(Operator& op,int pord, bool verbose) :
    ord_(pord),
    comm_(Communicator::instance()),
    op_(op),
    expl_(op),
    ode_(0)
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
  
  ~ExplTimeStepperBase() { delete ode_; }
  
protected:
  int ord_;
  Communicator & comm_;
  const Operator& op_;
  OperatorWrapper<Operator> expl_;
  DuneODE::ODESolver* ode_;
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
    BaseType(op,pord,verbose),
    timeProvider_(tp)
  {
    double cflLocal = (double (pord+1)/(pord+2));
    // maximal allowed cfl number 
    tp.provideCflEstimate(cflLocal); 
    assert( tp.cfl() <= 1.0 );
  }

  //! destructor 
  virtual ~ExplicitOdeSolver() {}
  
  //! solve system 
  void solve(DestinationType& U0) 
  {
    // if dt has not been set yet 
    if ( timeProvider_.notInitialized() ) 
    {
      DestinationType tmp("TMP",this->op_.space());
      this->op_(U0,tmp);
      return ;
    }
    
    // get dt 
    double dt = timeProvider_.deltaT();
    // should be larger then zero 
    assert( dt > 0.0 );
    
    // get time 
    double time = timeProvider_.time();

    // get leakPointer 
    double* u = U0.leakPointer();
    
    // call ode solver 
    const bool convergence = this->ode_->step(time, dt , u);

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ExplTimeStepper! \n";
      abort();
    }
  }

 private:
  const TimeProvider& timeProvider_;
};

template<class Operator>
class ExplTimeStepper : public TimeProvider, 
                        public ExplTimeStepperBase<Operator>  
{
  typedef ExplTimeStepperBase<Operator> BaseType;
  typedef typename Operator:: DestinationType ::
    DiscreteFunctionSpaceType :: GridType :: Traits ::
    CollectiveCommunication DuneCommunicatorType; 
 public:
  ExplTimeStepper(Operator& op,int pord, double cfl, bool verbose = false) :
    TimeProvider(0.0,cfl),
    BaseType(op,pord,verbose),
    tp_(this->op_.space().grid().comm(), *this ),
    savestep_(1),
    savetime_(0.0)
  {
    op.timeProvider(this);
    assert( this->cfl() <= 1.0 );
  }
  
  double solve(typename Operator::DestinationType& U0) 
  {
    if( tp_.notInitialized() ) 
    {
      typename Operator::DestinationType tmp("TMP",this->op_.space());
      this->op_(U0,tmp);
        
      // global min of dt 
      tp_.syncTimeStep(); 
    }
    
    tp_.resetTimeStepEstimate();
    
    // get dof array 
    double* u=U0.leakPointer();

    assert( tp_.deltaT() > 0.0 );
    
    const bool convergence = this->ode_->step(tp_.time(), tp_.deltaT() , u);

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ExplTimeStepper! \n";
      abort();
    }

    // global min of dt 
    tp_.syncTimeStep();
    
    // set new time to t + deltaT 
    tp_.augmentTime();

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
 public:
  ImplTimeStepperBase(Operator& op,int pord, bool verbose) :
    ord_(pord),
    comm_(Communicator::instance()),
    op_(op),
    impl_(op),
    ode_(0),
    linsolver_(comm_,cycle)
  {
    linsolver_.set_tolerance(1.0e-8);
    linsolver_.set_max_number_of_iterations(1000);
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
    ode_->set_tolerance(1.0e-8);
    if( verbose ) 
    {
      ode_->IterativeSolver::set_output(cout);
      ode_->DynamicalObject::set_output(cout);
    }
  }
  
  //! destructor 
  ~ImplTimeStepperBase() {delete ode_;}
  
protected:
  int ord_;
  Communicator & comm_;   
  const Operator& op_;
  OperatorWrapper<Operator> impl_;
  DuneODE::DIRK* ode_;
  DuneODE::GMRES linsolver_;
  enum { cycle = 20 };
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
    BaseType(op,pord,verbose),
    timeProvider_(tp)
  {
  }

  virtual ~ImplicitOdeSolver() {}
  
  void solve(DestinationType& U0) 
  {
    // for first call only calculate time step estimate 
    if( timeProvider_.notInitialized() ) 
    {
      DestinationType tmp("tmp",this->op_.space());
      this->op_(U0,tmp);
      return ;
    }

    double dt   = timeProvider_.deltaT();
    assert( dt > 0.0 );
    double time = timeProvider_.time();
    double* u = U0.leakPointer();
    const bool convergence = this->ode_->step(time , dt , u);

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ImplicitOdeSolver! \n";
      abort();
    }
  }

private:
  const TimeProvider& timeProvider_;
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
    BaseType(op,pord,verbose),
    tp_(this->op_.space().grid().comm(),*this),
    savestep_(1),
    savetime_(0.0)
  {
    op.timeProvider(this);
  }
  
  double solve(DestinationType& U0) 
  {
    if ( tp_.notInitialized() )     
    {
      DestinationType tmp("tmp",this->op_.space());

      this->op_(U0,tmp);

      // calculate global min of dt 
      tp_.syncTimeStep();
    }

    // reset dt estimate
    tp_.resetTimeStepEstimate();
    double* u=U0.leakPointer();
    
    assert( tp_.deltaT() > 0.0 );
    
    const bool convergence =  this->ode_->step(tp_.time() , tp_.deltaT() , u);

    assert(convergence);
    if(!convergence) 
    {
      std::cerr << "No Convergence of ImplTimeStepper! \n";
      abort();
    }

    // calculate global min of dt 
    tp_.syncTimeStep();

    // calls setTime(t+ this->dt_);
    tp_.augmentTime(); 
    
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
    expl_(op_expl),
    impl_(op_impl),
    ode_(0),
    linsolver_(comm_,cycle),
    tp_(this->op_.space().grid().comm(),*this),
    savetime_(0.0), savestep_(1)
  {
    op_expl.timeProvider(this);
    linsolver_.set_tolerance(1.0e-8);
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
    ode_->set_tolerance(1.0e-8);
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

    if ( tp_.notInitialized() ) 
    {
      typename OperatorExpl::DestinationType tmp("tmp",opexpl_.space());
      opexpl_(U0,tmp);
      
      // calculate global min of dt 
      tp_.syncTimeStep();
    }
    
    tp_.resetTimeStepEstimate();

    double* u=U0.leakPointer();
    
    const bool convergence = ode_->step( tp_.time(), tp_.deltaT() , u);
    assert(convergence);

    // calculate global min of dt 
    tp_.syncTimeStep();

    // calls setTime(t+dt_);
    tp_.augmentTime();

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
  OperatorWrapper<OperatorImpl> impl_;
  OperatorWrapper<OperatorExpl> expl_;
  DuneODE::SIRK* ode_;
  DuneODE::GMRES linsolver_;
  enum { cycle = 20 };
  // TimeProvider with communicator 
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
};

template<class Operator>
class ExplRungeKutta : public TimeProvider 
{
 public:
  enum {maxord=10};
  typedef typename Operator::SpaceType SpaceType;
  typedef typename Operator::DestinationType DestinationType;
  typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
 private:
  std::vector< std::vector<double> > a;
  std::vector<double> b;
  std::vector<double> c;
  
  int ord_;
  std::vector<DestinationType*> Upd;
public:
  ExplRungeKutta(Operator& op,int pord,double cfl, bool verbose = true ) :
    TimeProvider(0.0,cfl),
    op_(op),
    tp_(this->op_.space().grid().comm(),*this),
    ord_(pord), Upd(0),
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
    assert(ord_>0);
    a.resize(ord_);
    for (int i=0;i<ord_;i++)
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
    for (int i=0;i<ord_;i++)
    {
      Upd.push_back(new DestinationType("URK",op_.space()) );
    }
    
    Upd.push_back(new DestinationType("Ustep",op_.space()) );
  }

  ~ExplRungeKutta()
  {
    for(size_t i=0; i<Upd.size(); ++i) 
      delete Upd[i];
  }
  
  double solve(typename Operator::DestinationType& U0) 
  {
    tp_.resetTimeStepEstimate();
    
    double t = tp_.time();
    // Compute Steps
    op_(U0, *(Upd[0]));
    
    // global min of dt 
    tp_.syncTimeStep();
    
    // get cfl * timeStepEstimate 
    double dt = tp_.deltaT();

    for (int i=1;i<ord_;i++) 
    {
      (Upd[ord_])->assign(U0);
      for (int j=0; j<i ; j++) 
      {
        (Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      }

      setTime( t + c[i]*dt );

      op_(*(Upd[ord_]),*(Upd[i]));
    }

    // Perform Update
    for (int j=0;j<ord_;j++) 
    {
      U0.addScaled(*(Upd[j]),(b[j]*dt));
    }
    
    // calls setTime ( t + dt ); 
    tp_.setTime( t + dt );
    
    return tp_.time();
  }
  
  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ExplRungeKutta, steps: " << ord_ << "\n\n";
    ofs << "                cfl: " << tp_.cfl() << "\\\\\n\n";
    ofs.close();
    op_.printmyInfo(filename);
  }

 private:
  const Operator& op_;
  // TimeProvider with communicator 
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  int savestep_;
  double savetime_;
};

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
