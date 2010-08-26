#ifndef RUNGEKUTTA_SOLVER_HH
#define RUNGEKUTTA_SOLVER_HH

// include all used headers before, that they do not appear in DuneODE

//- system includes 
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

//- Dune includes 
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/timeprovider.hh>

namespace DuneODE 
{

using namespace Dune;
using namespace std;

/** @addtogroup ODESolver
 
    A variety of runge kutta have so far been implemented based
    on the \c pardg package. 
    These include explicit ssp methods upto order 4
    through Dune::ExplicitRungeKuttaSolver and Dune::ExplicitTimeStepper
    (where the second version uses \c pardg).
    Furthermore implicit methods and IMEX schemes up to order 
    4 can be used through Dune::ImplicitOdeSolver
    and Dune::SemiImplicitOdeSolver. 
    Each of these classes require an operator defining the right hand
    side for its construction; the template argument is a
    Dune::DiscreteFunctionSpace and the operator has to be derived
    from the Dune::SpaceOperatorInterface<DiscreteFunctionSpace>.
    Further arguments in the constructor are an instance of the
    Dune::Timeprovider and the requested order of the scheme which can
    be chosen during runtime.

    The management of the simulation time is performed by an instance of the
    Dune::TimeProvider class. The time at which the space operator is
    to be evaluated is passed to the operator in the \c setTime method on
    the Dune::SpaceOperatorInterface and it is assumed that
    through the method \c timeStepEstimate on the
    Dune::SpaceOperatorInterface an upper estimate for the time
    step is provided satifying at the least the stability
    restriction of the forward euler method. This estimate is then 
    multiplied by an additional factor depending on the ode solver used
    and then passed to the instance of the Dune::TimeProvider.

    \remarks 
    The interface for ODE solvers is defined by the class
    OdeSolverInterface. The interface for discretization operators
    working with the OdeSolvers is described by the class SpaceOperatorInterface.
 @{
 **/

/** \brief Interface class for ODE Solver. */ 
template <class DestinationImp>
class OdeSolverInterface 
{
protected:
  //! constructor
  OdeSolverInterface () {}    
public:
  //! type of destination 
  typedef DestinationImp DestinationType;

  //! destructor 
  virtual ~OdeSolverInterface () {}
  
  /** \brief initialize solver 
      \param[in] arg argument to apply internal operator once for intial time step estimate
  */
  virtual void initialize(const DestinationType& arg) = 0;
  
  /** \brief solve \f$\partial_t u = L(u)\f$ where \f$L\f$ is the internal operator.   
      \param[in] u unknown to solve for 
  */
  virtual void solve(DestinationType& u)
  {
    std::cerr << "OdeSolverInterface::solve(DestinationType) should not be used." 
              << std::endl;
    abort();
  }

  virtual void solve(DestinationType& u, int& newton_iterations, int& ils_iterations,
                     int& max_newton_iterations, int& max_ils_iterations)
  {
    std::cerr << "OdeSolverInterface::solve(DestinationType,int&,int&) should not be used." 
              << std::endl;
    abort();
  }

  /** \brief print description of ODE solver to out stream */
  virtual void description(std::ostream&) const = 0;
};

/** \brief Exlicit RungeKutta ODE solver. */
template<class DestinationImp>
class ExplicitRungeKuttaSolver : 
  public OdeSolverInterface<DestinationImp>
{
public:
  typedef DestinationImp DestinationType; 
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
protected:
  std::vector< std::vector<double> > a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<DestinationType*> Upd;
  const int ord_;

public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] tp TimeProvider 
    \param[in] pord polynomial order 
    \param[in] verbose verbosity 
  */
  ExplicitRungeKuttaSolver(OperatorType& op, 
                           TimeProviderBase& tp, 
                           const int pord, 
                           bool verbose = true ) :
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
    
    switch (ord_)
    {
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

  //! destructor 
  ~ExplicitRungeKuttaSolver()
  {
    for(size_t i=0; i<Upd.size(); ++i) 
      delete Upd[i];
  }

  //! apply operator once to get dt estimate 
  void initialize(const DestinationType& U0)
  {
    if( ! initialized_ ) 
    {
      // Compute Steps
      op_(U0, *(Upd[0]));
      initialized_ = true;
      
      // provide operators time step estimate 
      tp_.provideTimeStepEstimate( op_.timeStepEstimate() );
    }
  }
  
  //! solve the system 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ExplicitRungeKuttaSolver wasn't initialized before first call!");
    }

    // get cfl * timeStepEstimate 
    const double dt = tp_.deltaT();
    // get time 
    const double t = tp_.time();

    // set new time 
    op_.setTime( t );

    // Compute Steps
    op_(U0, *(Upd[0]));
    
    // provide operators time step estimate 
    tp_.provideTimeStepEstimate( op_.timeStepEstimate() );

    for (int i=1; i<ord_; ++i) 
    {
      (Upd[ord_])->assign(U0);
      for (int j=0; j<i ; ++j) 
      {
        (Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      }

      // set new time 
      op_.setTime( t + c[i]*dt );

      // apply operator 
      op_( *(Upd[ord_]), *(Upd[i]) );
      
      // provide operators time step estimate 
      tp_.provideTimeStepEstimate( op_.timeStepEstimate() );
    }

    // Perform Update
    for (int j=0; j<ord_; ++j) 
    {
      U0.addScaled(*(Upd[j]),(b[j]*dt));
    }
  }

  void description(std::ostream& out) const
  {
    out << "ExplRungeKutta, steps: " << ord_
        //<< ", cfl: " << this->tp_.factor() 
        << "\\\\" <<std::endl;
  }

protected:
  // operator to solve for 
  OperatorType& op_;
  // time provider 
  TimeProviderBase& tp_;
  // init flag 
  bool initialized_;
};

/** @} **/
} // end namespace DuneODE 
#endif
