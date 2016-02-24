#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

// inlcude all used headers before, that they don not appear in DuneODE
#ifdef HUNDERT
#endif
//- system includes
#include <iostream>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cassert>
#include <sys/times.h>

//- dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/odesolverinterface.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/solver/rungekutta/timestepcontrol.hh>

// include headers of PARDG
#include "pardg.hh"

namespace DuneODE
{

  using namespace Dune;
  using namespace Fem;
  using namespace std;

#ifdef USE_PARDG_ODE_SOLVER
  // see rungekutta/timestepcontrol.hh for the adaptive control of the clf for implicit solvers
  struct ODEParameters : public ImplicitRungeKuttaSolverParameters
  {
    using ImplicitRungeKuttaSolverParameters :: keyPrefix_;
    using ImplicitRungeKuttaSolverParameters :: tolerance ;

    ODEParameters( const std::string keyPrefix = "fem.ode." )
      : ImplicitRungeKuttaSolverParameters( keyPrefix )
    {}

    // choice of linear solver for the implicit ODE solver
    virtual PARDG::IterativeLinearSolver *linearSolver(PARDG::Communicator & comm) const
    {
      static const std::string methodTypeTable[] = { "gmres", "cg", "bicgstab" };
      const int method = Parameter::getEnum( keyPrefix_ + "linearsolver", methodTypeTable, 0 );

      PARDG::IterativeLinearSolver *solver = nullptr;
      switch( method )
      {
      case 0:
        solver = new PARDG::GMRES( comm, Parameter::getValue< int >( keyPrefix_ + "gmres.cycles", 15 ) );
        break;

      case 1:
        solver = new PARDG::CG( comm );
        break;

      case 2:
        solver = new PARDG::BICGSTAB( comm );
        break;
      }
      if( !solver )
        DUNE_THROW( InvalidStateException, "Unable to create linear solver." );

      // tolerance for the linear solver
      const double defaulTol = tolerance() * 1e-2 ;
      double tol = Parameter::getValue< double >( keyPrefix_ + "solver.tolerance" , defaulTol );
      std::string key( keyPrefix_ + "solver.errormeasure" );
      PARDG::set_tolerance(*solver,tol, key.c_str() );
      // max iterations that the linear solver should do
      int maxIter = Parameter::getValue< int >( keyPrefix_ + "solver.iterations" , 1000 );
      solver->set_max_number_of_iterations(maxIter);
      return solver;
    }

    // overload clone method
    virtual ODEParameters* clone () const { return new ODEParameters( *this ); }
  };


  template <class Operator>
  class OperatorWrapper : public PARDG::Function
  {
    // type of discrete function
    typedef typename Operator::DestinationType DestinationType;
   public:
    //! constructor
    OperatorWrapper(Operator& op)
      : op_(op)
    {}

    //! apply operator applies space operator and creates temporary
    //! discrete function using the memory from outside
    void operator()(const double *u, double *f, int i = 0)
    {
      // set actual time of iteration step
      op_.setTime( this->time() );

      // call operator apply
      op_( u, f );
    }

    //! return size of argument
    int dim_of_argument(int i = 0) const
    {
      assert( i == 0 );
      return op_.size();
    }

    //! return size of destination
    int dim_of_value(int i = 0) const
    {
      assert( i == 0 );
      return op_.size();
    }

    //! return reference to real operator
    const Operator& op() const { return op_; }

  protected:
    // operator to call
    Operator& op_;
  };


  template <class Operator>
  class LimiterWrapper : public PARDG::Limiter
  {
    // type of discrete function
    typedef typename Operator::DestinationType DestinationType;
   public:
    //! constructor
    LimiterWrapper(Operator& op)
      : op_(op)
    {}

    void operator()( double* ) { abort(); }

    //! apply operator applies space operator and creates temporary
    //! discrete function using the memory from outside
    void operator()(const double *u, double *f)
    {
      // call operator apply
      op_.limit( u, f );
    }

    //! return reference to real operator
    const Operator& op() const { return op_; }

  protected:
    // operator to call
    Operator& op_;
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
    typedef PARDGSpaceOperatorInterface<DestinationType> OperatorType;

    //! type of monitor class
    typedef typename OdeSolverInterface<DestinationImp> :: MonitorType  MonitorType;
  protected:
    //! constructor
    ParDGOdeSolverBase( TimeProviderBase& tp,
                        const int order ) :
      timeProvider_( tp ),
      comm_(PARDG::Communicator::instance()),
      order_( order ),
      initialized_( false ),
      odeSolver_( 0 )
    {
    }

    void initializeOdeSolver ()
    {
      if( !initialized_ )
      {
        if( !odeSolver_ )
          odeSolver_ = createOdeSolver();
        initialized_ = true;
      }
    }

  public:
    // initialize time step size
    void initialize ( const DestinationType& U0 )
    {
      // initialized dt on first call
      if ( ! initialized_ )
      {
        // create instance of ode solver
        initializeOdeSolver();

        // apply operator once
        spaceOperator().initializeTimeStepSize( U0 );

        // set initial time step
        timeProvider_.provideTimeStepEstimate( spaceOperator().timeStepEstimate() );
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

  public:
    void description(std::ostream& out) const
    {
      out << name() << ", steps: " << order_;
      // would be nice to add info on cfl number, but TimeProviderBase doesn't have this
          //<< ", cfl: " << this->timeProvider_.cfl() << "\\\\" <<std::endl;
    }

  protected:
    //! return name of ode solver
    virtual std::string name()  const = 0;

    //! return ode solver object
    virtual PARDG::ODESolver* createOdeSolver() = 0;

    //! return reference to discretization operator for time step initialization
    virtual const OperatorType& spaceOperator() const = 0;

  protected:
    TimeProviderBase& timeProvider_;
    PARDG::Communicator & comm_;
    const int order_;
    bool initialized_;

  private:
    // use method odeSolver for access
    PARDG::ODESolver* odeSolver_;
  };

  ///////////////////////////////////////////////////////
  //
  //  --ExplicitOdeSolver
  //
  ///////////////////////////////////////////////////////
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

  public:
    using BaseType :: initialize;
    using BaseType :: solve ;

    typedef typename BaseType :: OperatorType    OperatorType;
    typedef typename BaseType :: DestinationType DestinationType;
    typedef typename BaseType :: MonitorType     MonitorType;

    //! constructor
    ExplicitOdeSolver(OperatorType& op,
                      TimeProviderBase &tp,
                      const int order, bool verbose = false) :
      BaseType(tp, order ),
      expl_( op ),
      verbose_( verbose )
    {}

    //! destructor
    virtual ~ExplicitOdeSolver() {}

    //! solve system
    void solve( DestinationType& U0,
                MonitorType& monitor )
    {
      // no nonlinear system to solve for time step update
      monitor.reset() ;

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
      const bool convergence = odeSolver().step(time, dt, u);

      // set time step estimate of operator
      timeProvider_.provideTimeStepEstimate( expl_.op().timeStepEstimate() );

      assert( convergence );
      if(! convergence )
      {
        timeProvider_.invalidateTimeStep();
        if( comm_.id() == 0 )
          std::cerr << "No Convergence of ExplicitOdeSolver! \n";
      }
    }

  protected:
    //! return name of ode solver
    virtual std::string name() const { return "ExplicitOdeSolver"; }

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
                  if( comm_.id() == 0 )
                  {
                    std::cerr << "Runge-Kutta method of this order not implemented.\n"
                              << "Using 7-stage Bulirsch-Stoer scheme.\n"
                              << std::endl;
                  }
      }

      // only set output when general verbose mode is enabled
      // (basically to avoid output on every rank)
      if( verbose_ && Parameter :: verbose() )
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

  public:
    using BaseType :: initialize;
    using BaseType :: solve ;

    typedef typename BaseType :: OperatorType    OperatorType;
    typedef typename BaseType :: DestinationType DestinationType;
    typedef typename BaseType :: MonitorType     MonitorType;

    ImplicitOdeSolver(OperatorType& op,
                      TimeProviderBase& tp,
                      const int order,
                      const ODEParameters& parameter = ODEParameters() ) :
      BaseType(tp, order),
      impl_( op ),
      linsolver_( 0 ),
      param_( parameter.clone() ),
      verbose_( parameter.verbose() ),
      cfl_( parameter.cflStart() ),
      cflMax_( parameter.cflMax() )
    {
    }

    virtual ~ImplicitOdeSolver()
    {
      delete linsolver_; linsolver_ = 0;
      delete param_;     param_ = 0;
    }

  protected:
    //! return name of ode solver
    virtual std::string name() const { return "ImplicitOdeSolver"; }

    //! return reference to parameter class
    const ODEParameters& parameter() const { assert( param_ ); return *param_; }

    //! for initialization
    virtual const OperatorType& spaceOperator() const { return impl_.op(); }

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
        case 4: odeSolver = new PARDG::DIRK34 (comm_, impl_); break;
        default : if( comm_.id() == 0 )
                  {
                    std::cerr << "Runge-Kutta method of this order not implemented"
                              << std::endl;
                  }
                  abort();
      }

      assert( odeSolver );

      odeSolver->set_linear_solver( *linsolver_ );
      odeSolver->set_tolerance( parameter().tolerance() );
      odeSolver->set_max_number_of_iterations( parameter().iterations() );

      // only set output when general verbose mode is enabled
      // (basically to avoid output on every rank)
      if( verbose_ == ODEParameters :: fullVerbosity &&
          Parameter :: verbose() )
      {
        odeSolver->IterativeSolver::set_output(cout);
        odeSolver->DynamicalObject::set_output(cout);
        linsolver_->IterativeSolver::set_output(cout);
      }
      return odeSolver;
    }

    virtual int numberOfIterations()
    {
      return static_cast<PARDG::DIRK&> (odeSolver()).number_of_iterations();
    }

  public:
    //! solve
    void solve( DestinationType& U0, MonitorType& monitor )
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

      const bool convergence =
        odeSolver().step( time, dt, u,
                          monitor.newtonIterations_,
                          monitor.linearSolverIterations_,
                          monitor.maxNewtonIterations_,
                          monitor.maxLinearSolverIterations_ );

      assert( linsolver_->number_of_iterations() == monitor.linearSolverIterations_ );

      double factor( 1 );
      bool changed =
        parameter().cflFactor( impl_.op().timeStepEstimate(), spaceOperator().timeStepEstimate(),
                               monitor.linearSolverIterations_, convergence, factor );

      if( (factor >= std::numeric_limits< double >::min()) &&
          (factor <= std::numeric_limits< double >::max()) )
      {
        // only apply when factor is small or max cfl was not reached yet
        if( factor < 1.0 || cfl_ <= cflMax_ )
        {
          cfl_ *= factor;
        }
        else
          changed = false ;
      }
      else
        DUNE_THROW( InvalidStateException, "invalid cfl factor: " << factor );

      if (convergence)
      {
        // timeProvider_.provideTimeStepEstimate( cfl_ * spaceOperator().timeStepEstimate() );
        timeProvider_.provideTimeStepEstimate( timeStepEstimate(cfl_) );

        if( changed && (verbose_ >= ODEParameters :: cflVerbosity ) && (MPIManager::rank() == 0) )
        {
          derr << " New cfl number is: "<< cfl_ << ", iterations per time step("
               << "ILS: " << monitor.linearSolverIterations_
               << ", Newton: " << monitor.newtonIterations_
               << ")"
               << std::endl;
        }
      }
      else
      {
        if( factor >= 1.0 )
          DUNE_THROW( InvalidStateException, "Invalid cfl factor (no convergence): " << factor );
        timeProvider_.provideTimeStepEstimate( factor * dt );
        timeProvider_.invalidateTimeStep();

        // output only in verbose mode
        if( (verbose_ >= ODEParameters :: noConvergenceVerbosity) && (MPIManager::rank() == 0) )
          derr << "No convergence: New cfl number is " << cfl_ << std::endl;
      }

      linsolver_->reset_number_of_iterations();
    }

  protected:
    virtual double timeStepEstimate(double cfl)
    {
      return cfl*spaceOperator().timeStepEstimate();
    }

    OperatorWrapper<OperatorType> impl_;
    PARDG::IterativeLinearSolver* linsolver_;
    const ODEParameters* param_;
    const int verbose_;
    double cfl_;
    const double cflMax_;
  }; // end ImplicitOdeSolver


  ///////////////////////////////////////////////////////
  //
  //  --SemiImplicitOdeSolver
  //
  ///////////////////////////////////////////////////////
  template<class DestinationImp>
  class SemiImplicitOdeSolver :
    public ImplicitOdeSolver< DestinationImp >
  {
    typedef ImplicitOdeSolver< DestinationImp > BaseType;
  protected:
    using BaseType :: order_;
    using BaseType :: comm_;
    using BaseType :: initialized_;
    using BaseType :: timeProvider_;
    using BaseType :: odeSolver ;
    using BaseType :: linsolver_ ;
    using BaseType :: parameter ;
    using BaseType :: impl_ ;
    using BaseType :: verbose_ ;
    using BaseType :: cfl_;

  public:
    typedef typename BaseType :: OperatorType    OperatorType;
    typedef typename BaseType :: DestinationType DestinationType;

    SemiImplicitOdeSolver(OperatorType& explOp,
                          OperatorType& implOp,
                          TimeProviderBase& tp,
                          const int order,
                          const ODEParameters& parameter = ODEParameters() ) :
      BaseType( implOp, tp, order, parameter ),
      expl_( explOp ),
      limiter_( explOp )
    {
    }

    // initialize time step size
    void initialize ( const DestinationType& U0 )
    {
      // initialized dt on first call
      if( !initialized_ )
      {
        BaseType::initializeOdeSolver();

        // apply operators once
        expl_.op().initializeTimeStepSize( U0 );
        impl_.op().initializeTimeStepSize( U0 );

        // set initial time step
        double dtEst = std::numeric_limits< double >::max();
        parameter().initTimeStepEstimate( expl_.op().timeStepEstimate(), impl_.op().timeStepEstimate(), dtEst, cfl_ );
        timeProvider_.provideTimeStepEstimate( dtEst );
      }
    }


  protected:
    //! return name of ode solver
    virtual std::string name() const { return "SemiImplicitOdeSolver"; }

    //! for initialization
    const OperatorType& spaceOperator() const { return expl_.op(); }

    //! create implicit ode solver
    PARDG::ODESolver* createOdeSolver()
    {
      linsolver_ = parameter().linearSolver( comm_ );
      assert( linsolver_ );

      PARDG::SIRK* odeSolver = 0;
      switch ( order_ )
      {
        case 1: odeSolver = new PARDG::SemiImplicitEuler(comm_, impl_, expl_); break;
        case 2: odeSolver = new PARDG::IMEX_SSP222      (comm_, impl_, expl_); break;
        case 3: odeSolver = new PARDG::SIRK33           (comm_, impl_, expl_); break;
                //odeSolver = new PARDG::IMEX_ARK34       (comm_, impl_, expl_); break;
        case 4: odeSolver = new PARDG::IERK45           (comm_, impl_, expl_); break;
                //odeSolver = new PARDG::IMEX_ARK46       (comm_, impl_, expl_); break;
        default : if( comm_.id() == 0 )
                  {
                    std::cerr << "Runge-Kutta method of this order not implemented"
                              << std::endl;
                  }
                  abort();
      }

      assert( odeSolver );

      odeSolver->set_linear_solver(*linsolver_);
      odeSolver->set_tolerance( parameter().tolerance() );
      odeSolver->set_max_number_of_iterations( parameter().iterations() );

      if( expl_.op().hasLimiter() )
        odeSolver->set_expl_limiter( limiter_ );

      // only set output when general verbose mode is enabled
      // (basically to avoid output on every rank)
      if( verbose_ == ODEParameters :: fullVerbosity &&
          Parameter :: verbose() )
      {
        odeSolver->IterativeSolver::set_output(cout);
        odeSolver->DynamicalObject::set_output(cout);
        // linsolver_->set_output(cout);
      }
      return odeSolver;
    }

    virtual int numberOfIterations()
    {
      return static_cast<PARDG::SIRK&> (odeSolver()).number_of_iterations();
    }

  protected:
    virtual double timeStepEstimate(double cfl)
    {
      // take the minimum of the explicit part and the implicit part scaled by the cfl number
      // return spaceOperator().timeStepEstimate();
      return std::min( spaceOperator().timeStepEstimate(),
                       cfl * impl_.op().timeStepEstimate() );
    }
    OperatorWrapper<OperatorType> expl_;
    LimiterWrapper <OperatorType> limiter_;

  }; // end SemiImplicitOdeSolver

#endif // USE_PARDG_ODE_SOLVER

  /**
   @}
  **/

} // namespace DuneODE

#undef USE_EXTERNAL_BLAS
#endif // #ifndef RUNGEKUTTA_ODE_SOLVER_HH
