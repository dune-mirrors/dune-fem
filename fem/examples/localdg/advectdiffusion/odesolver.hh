#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

#include <iostream>
#include <cmath>
#include "../../../misc/timeutility.hh"

 #include "ode/function.hpp"
 #include "ode/ode_solver.hpp"
 #include "ode/linear_solver.hpp"
 #include "ode/bulirsch_stoer.cpp"  
 #include "ode/iterative_solver.cpp"  
 #include "ode/ode_solver.cpp"     
 #include "ode/sirk.cpp"   
 // #include "ode/communicator.cpp"    
 #include "ode/matrix.cpp"            
 #include "ode/qr_solver.cpp"   
 #include "ode/ssp.cpp"
 #include "ode/dirk.cpp"     
 #include "ode/runge_kutta.cpp"
 #include "ode/vector.cpp"


namespace DuneODE {
  using namespace Dune;
  using namespace std;

template <class Operator>
 class OperatorWrapper : public DuneODE::Function {
 public:
  OperatorWrapper(const Operator& op) : op_(op) {}
  void operator()(const double *u, double *f, int i = 0) {
    typename Operator::DestinationType arg("ARG",op_.space(),u);
    typename Operator::DestinationType dest("DEST",op_.space(),f);
    // op_.setTime(time());
    op_(arg,dest);
  }
  int dim_of_argument(int i = 0) const 
  { 
    if (i==0) return op_.space().size();
    else assert(0);
  }
  int dim_of_value(int i = 0) const 
  { 
    if (i==0) return op_.space().size();
    else assert(0);
  }
 private:
  const Operator& op_;
};
template<class Operator>
class ExplTimeStepper : public TimeProvider {
 public:
  ExplTimeStepper(Operator& op,int pord,double cfl) :
    op_(op),
    expl_(op),
    ode_(0),
    cfl_(cfl),
    dt_(-1.0),
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
    switch (pord) {
    case 1: ode_=new ExplicitEuler(comm,expl_); break;
    case 2: ode_=new ExplicitTVD2(comm,expl_); break;
    case 3: ode_=new ExplicitTVD3(comm,expl_); break;
    case 4: ode_=new ExplicitRK4(comm,expl_); break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
			<< std::endl;
      abort();
    }
    ode_->DynamicalObject::set_output(cout);
  }
  ~ExplTimeStepper() {delete ode_;}
  double solve(typename Operator::DestinationType& U0) {
    if (dt_<0) {
      typename Operator::DestinationType tmp("TMP",op_.space());
      op_(U0,tmp);
      dt_=cfl_*timeStepEstimate();
    }
    resetTimeStepEstimate();
    double t=time();
    double* u=U0.leakPointer();
    const bool convergence = ode_->step(t, dt_, u);
    assert(convergence);
    setTime(t+dt_);
    dt_=cfl_*timeStepEstimate();
    return time();
  }
  void printGrid(int nr, 
		 const typename Operator::DestinationType& U) {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
 private:
  DuneODE::Communicator comm;
  const Operator& op_;
  OperatorWrapper<Operator> expl_;
  DuneODE::ODESolver* ode_;
  double cfl_;
  double dt_;
  int savestep_;
  double savetime_;
};
template<class Operator>
class ImplTimeStepper : public TimeProvider {
 public:
  ImplTimeStepper(Operator& op,int pord,double cfl) :
    op_(op),
    impl_(op),
    ode_(0),
    linsolver_(comm),
    cfl_(cfl),
    dt_(-1.0),
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
    linsolver_.set_tolerance(1.0e-8);
    linsolver_.set_max_number_of_iterations(1000);
    switch (pord) {
    case 1: ode_=new ImplicitEuler(comm,impl_); break;
    case 2: ode_=new Gauss2(comm,impl_); break;
    case 3: ode_=new DIRK3(comm,impl_); break;
      //case 4: ode_=new ExplicitRK4(comm,expl_); break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
			<< std::endl;
      abort();
    }
    ode_->set_linear_solver(linsolver_);
    ode_->set_tolerance(1.0e-8);
    ode_->IterativeSolver::set_output(cout);
    ode_->DynamicalObject::set_output(cout);
  }
  ~ImplTimeStepper() {delete ode_;}
  double solve(typename Operator::DestinationType& U0) {
    if (dt_<0) {
      typename Operator::DestinationType tmp("tmp",op_.space());
      op_(U0,tmp);
      dt_ = cfl_*timeStepEstimate();
    }
    resetTimeStepEstimate();
    double t=time();
    double* u=U0.leakPointer();
    const bool convergence = ode_->step(t, dt_, u);
    assert(convergence);
    setTime(t+dt_);
    dt_ = cfl_*timeStepEstimate();
    return time();
  }
  void printGrid(int nr, 
		 const typename Operator::DestinationType& U) {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
 private:
  DuneODE::Communicator comm;	  
  const Operator& op_;
  OperatorWrapper<Operator> impl_;
  DuneODE::DIRK* ode_;
  DuneODE::GMRES<20> linsolver_;
  double cfl_;
  double dt_;
  int savestep_;
  double savetime_;
};
template<class OperatorExpl,class OperatorImpl>
class SemiImplTimeStepper : public TimeProvider {
  typedef OperatorExpl Operator;
 public:
  SemiImplTimeStepper(OperatorExpl& op_expl,OperatorImpl& op_impl,
		      int pord,double cfl) :
    opexpl_(op_expl),
    opimpl_(op_impl),
    expl_(op_expl),
    impl_(op_impl),
    ode_(0),
    linsolver_(comm),
    cfl_(cfl),
    dt_(-1.0),
    savetime_(0.0), savestep_(1)
  {
    op_expl.timeProvider(this);
    linsolver_.set_tolerance(1.0e-8);
    linsolver_.set_max_number_of_iterations(1000);
    switch (pord) {
    case 1: ode_=new SemiImplicitEuler(comm,impl_,expl_); break;
    case 2: ode_=new IMEX_SSP222(comm,impl_,expl_); break;
    case 3: ode_=new SIRK33(comm,impl_,expl_); break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
			<< std::endl;
      abort();
    }
    ode_->set_linear_solver(linsolver_);
    ode_->set_tolerance(1.0e-8);
    ode_->IterativeSolver::set_output(cout);
    ode_->DynamicalObject::set_output(cout);
  }
  ~SemiImplTimeStepper() {delete ode_;}
  double solve(typename Operator::DestinationType& U0) {
    if (dt_<0) {
      typename OperatorExpl::DestinationType tmp("tmp",opexpl_.space());
      opexpl_(U0,tmp);
      dt_ = cfl_*timeStepEstimate();
    }
    resetTimeStepEstimate();
    double t=time();
    double* u=U0.leakPointer();
    const bool convergence = ode_->step(t, dt_, u);
    assert(convergence);
    setTime(t+dt_);
    dt_ = cfl_*timeStepEstimate();
    return time();
  }
  void printGrid(int nr, 
		 const typename Operator::DestinationType& U) {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,opexpl_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
 private:
  DuneODE::Communicator comm;	  
  const OperatorExpl& opexpl_;
  const OperatorImpl& opimpl_;
  OperatorWrapper<OperatorImpl> impl_;
  OperatorWrapper<OperatorExpl> expl_;
  DuneODE::SIRK* ode_;
  DuneODE::GMRES<20> linsolver_;
  double cfl_;
  double dt_;
  int savestep_;
  double savetime_;
};
template<class Operator>
class ExplRungeKutta : public TimeProvider {
 public:
  enum {maxord=10};
  typedef typename Operator::SpaceType SpaceType;
  typedef typename Operator::DestinationType DestinationType;
 private:
  double cfl_;
  double **a;
  double *b;
  double *c;
  int ord_;
  std::vector<DestinationType*> Upd;
public:
  ExplRungeKutta(Operator& op,int pord,double cfl) :
    op_(op),
    cfl_(cfl), ord_(pord), Upd(0),
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
    assert(ord_>0);
    a=new double*[ord_];
    for (int i=0;i<ord_;i++)
      a[i]=new double[ord_];
    b=new double [ord_];
    c=new double [ord_];
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
      Upd.push_back(new DestinationType("URK",op_.space()) );
    Upd.push_back(new DestinationType("Ustep",op_.space()) );
  }
  double solve(typename Operator::DestinationType& U0) {
    resetTimeStepEstimate();
    double t=time();
    // Compute Steps
    op_(U0,*(Upd[0]));
    double dt=cfl_*timeStepEstimate();
		cout << t << " " << dt << endl;
    for (int i=1;i<ord_;i++) {
      (Upd[ord_])->assign(U0);
      for (int j=0;j<i;j++) 
	(Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      setTime(t+c[i]*dt);
      op_(*(Upd[ord_]),*(Upd[i]));
      double ldt=cfl_*timeStepEstimate();
    }
    // Perform Update
    for (int j=0;j<ord_;j++) {
      U0.addScaled(*(Upd[j]),(b[j]*dt));
    }
    setTime(t+dt);
    return time();
  }
  void printGrid(int nr, 
		 const typename Operator::DestinationType& U) {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
 private:
  const Operator& op_;
  int savestep_;
  double savetime_;
};

}
#endif
