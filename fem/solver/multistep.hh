#ifndef MULTISTEP_SOLVER_HH
#define MULTISTEP_SOLVER_HH

// inlcude all used headers before, that they don not appear in DuneODE 

//- system includes 
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

//- Dune includes 
#include <dune/fem/misc/timeprovider.hh>
#include <dune/fem/solver/rungekutta.hh>
namespace DuneODE 
{

using namespace Dune;
using namespace std;

/** @ingroup ODESolver

 @{
 **/


/** \brief Base class for explicit multistep method. */
template<class Operator>
class ExplMultiStepBase 
{
public:
  typedef typename Operator::DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  // typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
private:
  std::vector< std::vector<double> > a;
  std::vector<double> b;
  std::vector<double> c;
  size_t steps_;
  double gamma_;
  std::vector<DestinationType*> Uj;
  std::vector<DestinationType*> Fj;
  std::vector<double> deltat_;
  bool msInit;
protected:  
  const int ord_;

public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] tp TimeProvider 
    \param[in] pord polynomial order
    \param[in] verbose verbosity 
  */
  ExplMultiStepBase(Operator& op, TimeProvider& tp, 
                    int pord, bool verbose = true ) :
    a(0),b(0),c(0)
    , steps_(3)
    , gamma_(1./3.)
    , Uj(0)
    , Fj(0)
    , deltat_(steps_)
    , msInit(false)
    , ord_(pord)
    , op_(op)
    , tp_(tp)
    , initialized_(false)
  {
    assert(ord_==2);
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
    // create memory for intermediate stages
    for (size_t i=0; i<steps_; ++i)
    {
      Fj.push_back(new DestinationType("FMS",op_.space()) );
    }
    Uj.push_back(new DestinationType("UMS",op_.space()) );
  }

  //! destructor 
  ~ExplMultiStepBase()
  {
    for(size_t i=0; i<Fj.size(); ++i) {
      delete Fj[i];
    }
    for(size_t i=0; i<Uj.size(); ++i) {
      delete Uj[i];
    }
  }

  //! apply operator once to get dt estimate 
  void initialize(const DestinationType& U0)
  {
    if( ! initialized_ ) 
    {
      // Compute Steps
      op_(U0, (*Fj[0]));
      initialized_ = true;
    }
  }
  
  //! solve the system 
  void solve(DestinationType& U0) 
  {
    if (!msInit) {
      // Startup phase using RungeKutta
      Uj[Uj.size()-1]->assign(U0);
      solveRK(U0);
      deltat_[Uj.size()-1] = tp_.deltaT();
      Uj.push_back(new DestinationType("UMS",op_.space()) );
      if (Uj.size()==steps_) {
        for (size_t i=0;i<steps_-1;i++) {
          op_(*(Uj[i]),(*Fj[i])); 
        }
        msInit=true;
      }
      std::cerr << "Startup with deltat = "
                << deltat_[0] << std::endl;
      return;
    }
    // now multistep

    // time might change 
    tp_.unlock();
    
    // get cfl * timeStepEstimate 
    deltat_[steps_-1] = tp_.deltaT();
    // get time 

    // Compute Steps
    Uj[steps_-1]->assign(U0);
    op_(*(Uj[steps_-1]), *(Fj[steps_-1]));

    // Perform Update
    double alpha[steps_];
    double beta[steps_];
    double w = deltat_[steps_-1]/deltat_[steps_-2];
    switch (steps_) {
    case 2: {
      double alpha2 = (1.+2.*gamma_*w)/(1.+w);
      alpha[1] = -((1.-2.*gamma_)*w-1.)/alpha2;
      alpha[0] = -((2.*gamma_-1.)*w*w/(1.+w))/alpha2;
      beta[1]  = (1.+gamma_*w)/alpha2*deltat_[steps_-1];
      beta[0]  = -(gamma_*w)/alpha2*deltat_[steps_-1];
    } break;
    case 3: {
      /* Shu
      alpha[0] = 1./4.;
      alpha[1] = 0.;
      alpha[2] = 3./4.;
      beta[0]  = 0.*deltat_[steps_-1];
      beta[1]  = 0.*deltat_[steps_-1];
      beta[2]  = 3./2.*deltat_[steps_-1];
      */
      double w0 = deltat_[1]/deltat_[0];
      double w1 = deltat_[2]/deltat_[1];
      double g = w/(w0+w1);
      // CFL < 1-g -> TVD
      // therefore need g<1!
      assert(g<1);
      alpha[0] = g*g;
      alpha[1] = 0.;
      alpha[2] = 1.-g*g; // > 0 due to CFL
      beta[0]  = 0.*deltat_[steps_-1];
      beta[1]  = 0.*deltat_[steps_-1];
      beta[2]  = (1+g)*deltat_[steps_-1]; // always > 0
    } break;
    case 4: {
      alpha[0] = 1./9.;
      alpha[1] = 0.;
      alpha[2] = 0.;
      alpha[3] = 8./9.;
      beta[0]  = 0.*deltat_[steps_-1];
      beta[1]  = 0.*deltat_[steps_-1];
      beta[2]  = 0.*deltat_[steps_-1];
      beta[3]  = 4./3.*deltat_[steps_-1];
    } break;
    }
    // Update
    U0 *= alpha[steps_-1];
    U0.addScaled(*(Fj[steps_-1]),beta[steps_-1]);
    for (size_t i=0;i<steps_-1;i++) {
      U0.addScaled(*(Uj[i]), alpha[i]);
      U0.addScaled(*(Fj[i]), beta[i]);
    }

    DestinationType* Utmp = Uj[0];
    DestinationType* Ftmp = Fj[0];
    for (size_t i=1;i<steps_;i++) {
      deltat_[i-1] = deltat_[i];
      Uj[i-1] = Uj[i];
      Fj[i-1] = Fj[i];
    }
    Uj[steps_-1] = Utmp;
    Fj[steps_-1] = Ftmp;
    // restore global time 
    tp_.lock();
  }
  //! solve the system 
  void solveRK(DestinationType& U0) 
  {
    DestinationType Uval("Utmp",op_.space());
    // time might change 
    tp_.unlock();
    
    // get cfl * timeStepEstimate 
    const double dt = tp_.deltaT();
    // get time 
    const double t = tp_.time();

    // Compute Steps
    op_(U0, *(Fj[0]));
    
    for (int i=1; i<ord_; ++i) 
    {
      Uval.assign(U0);
      for (int j=0; j<i ; ++j) 
      {
        Uval.addScaled(*(Fj[j]),(a[i][j]*dt));
      }

      // set new time 
      tp_.setTime( t + c[i]*dt );

      // apply operator 
      op_(Uval,*(Fj[i]));
    }

    // Perform Update
    for (int j=0; j<ord_; ++j) 
    {
      U0.addScaled(*(Fj[j]),(b[j]*dt));
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

/** \brief Exlicit multi step ODE solver that also behaves like a time
    stepper. */
template<class Operator>
class ExplMultiStep : public TimeProvider , 
                      public ExplMultiStepBase<Operator> 
{
  typedef ExplMultiStepBase<Operator> BaseType;
public:
  typedef typename Operator :: DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 

public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] pord polynomial order 
    \param[in] cfl cfl number 
    \param[in] verbose verbosity 
  */
  ExplMultiStep (Operator& op,int pord,double cfl, bool verbose = true ) :
    TimeProvider(0.0,cfl/double(pord)),
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
    ofs << "ExplMultiStep, steps: " << this->ord_ << "\n\n";
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

/** \brief Exlicit multi step ODE solver. */
template<class DestinationImp>
class ExplicitMultiStepSolver : 
  public OdeSolverInterface<DestinationImp> ,
  public ExplMultiStepBase<SpaceOperatorInterface<DestinationImp> >  
{
  typedef DestinationImp DestinationType; 
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef ExplMultiStepBase<OperatorType> BaseType;
 public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] tp TimeProvider 
    \param[in] pord polynomial order 
    \param[in] verbose verbosity 
  */
  ExplicitMultiStepSolver(OperatorType& op, TimeProvider& tp, int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {
    // CFL upper estimate 
    double cfl = 0.45 / (2.0 * pord+1) / double(pord);

    // maximal allowed cfl number 
    tp.provideCflEstimate(cfl); 
    assert( tp.cfl() <= 1.0 );

    if(verbose) 
    {
      std::cout << "ExplicitMultiStepSolver: cfl = " << tp.cfl() << "!\n";
    } 
  }

  //! destructor 
  virtual ~ExplicitMultiStepSolver() {}
  
  //! apply operator once to get dt estimate 
  void initialize(const DestinationType& U0)
  {
    BaseType :: initialize(U0);
  }
  
  //! solve system 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! this->initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ExplicitMultiStepSolver wasn't initialized before first call!");
    }

    // solve ode 
    BaseType :: solve(U0);
  }

private:
  TimeProvider& timeProvider_;
};

/** @} **/
} // end namespace DuneODE 
#endif
