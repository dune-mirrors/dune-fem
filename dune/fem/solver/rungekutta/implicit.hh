#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH

//- system includes 
#include <cassert>
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/rungekutta/basicimplicit.hh>
#include <dune/fem/solver/rungekutta/butchertable.hh>
#include <dune/fem/solver/timeprovider.hh>

namespace DuneODE 
{

  // ImplicitRungeKuttaSolverParameters
  // ----------------------------------

  struct ImplicitRungeKuttaSolverParameters
  : public Dune::Fem::LocalParameter< ImplicitRungeKuttaSolverParameters, ImplicitRungeKuttaSolverParameters >
  { 
    enum { noVerbosity = 0,  noConvergenceVerbosity = 1,
           cflVerbosity = 2, fullVerbosity = 3 };

    ImplicitRungeKuttaSolverParameters ()
    : min_it( Dune::Fem::Parameter::getValue< int >( "fem.ode.miniterations" , 14 ) ),
      max_it( Dune::Fem::Parameter::getValue< int >( "fem.ode.maxiterations" , 16 ) ),
      sigma( Dune::Fem::Parameter::getValue< double >( "fem.ode.cflincrease" , 1.1 ) )
    {}

    /** \brief tolerance for the non-linear solver (should be larger than the tolerance for
               the linear solver */
    virtual double tolerance () const
    {
      return Dune::Fem::Parameter::getValue< double >( "fem.ode.tolerance" , 1e-6 );
    }

    virtual int iterations() const
    {
      return Dune::Fem::Parameter::getValue< int >( "fem.ode.iterations" , 1000 );
    }

    /** \brief verbosity level ( none, noconv, cfl, full )  */
    virtual int verbose () const
    {
      static const std::string verboseTypeTable[] = { "none", "noconv", "cfl", "full" };
      return Dune::Fem::Parameter::getEnum( "fem.ode.verbose", verboseTypeTable, 0 );
    }

    virtual double cflStart () const
    {
      return Dune::Fem::Parameter::getValue< double >( "fem.ode.cflStart", 1 );
    }

    virtual double cflMax () const
    {
      return Dune::Fem::Parameter::getValue< double >( "fem.ode.cflMax" , std::numeric_limits< double >::max() );
    }

    double initialDeltaT ( double dt ) const
    {
      return std::min( Dune::Fem::Parameter::getValue< double >( "fem.ode.initialdt", 987654321 ), dt );
    }

    /** \brief return multiplication factor for the current cfl number
     *  \param[in] imOpTimeStepEstimate Time step estimate of the first ode solver
     *  \param[in] exOpTimeStepEstimate Time step estimate of the second ode solver
     *  \param[in] solver Iterative linear solver (ILS)
     *  \param[in] converged Convergence of the ILS
     *  \param[out] factor Multiplication factor for the current cfl number
     *
     *  \note Do not increase the cfl number of the implicit solver if its time step 
     *    estimate is already larger than the one of the explicit solver
     */
    virtual bool cflFactor( const double imOpTimeStepEstimate,
                            const double exOpTimeStepEstimate,
                            const int numberOfLinearIterations,
                            bool converged,
                            double &factor) const
    {
      const int iter = numberOfLinearIterations; 
      factor = 1.;
      bool changed = false;
      if (converged) 
      {
        if (iter < min_it) 
        {
          if( imOpTimeStepEstimate <= exOpTimeStepEstimate )
          {
            factor = sigma;
            changed = true;
          }
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

    virtual void initTimeStepEstimate ( const double dtEstExpl, const double dtEstImpl, double &dtEst, double &cfl ) const
    {
      // initial time step already set to explicit time step
      dtEst = dtEstExpl;

      // heuristics for initial CFL number
      cfl = 1.0;
      if( (dtEstImpl > 0) && (dtEstExpl > dtEstImpl) )
        cfl = dtEstExpl / dtEstImpl;
    }

    // number of minimal iterations that the linear solver should do 
    // if the number of iterations done is smaller then the cfl number is increased  
    int min_it;
    // number of maximal iterations that the linear solver should do 
    // if the number of iterations larger then the cfl number is decreased   
    int max_it;
    // factor for cfl number on increase (decrease is 0.5)
    double sigma;
  };



  // ImplicitRungeKuttaTimeStepControl
  // ---------------------------------

  class ImplicitRungeKuttaTimeStepControl
  {
    typedef ImplicitRungeKuttaTimeStepControl ThisType;

  public:
    typedef Dune::Fem::TimeProviderBase TimeProviderType;
    typedef ImplicitRungeKuttaSolverParameters ParametersType;

    explicit ImplicitRungeKuttaTimeStepControl ( TimeProviderType &timeProvider,
                                                 const ParametersType &parameters = ParametersType() )
    : timeProvider_( timeProvider ),
      parameters_( parameters.clone() ),
      cfl_( parameters_->cflStart() ),
      cflMax_( parameters_->cflMax() ),
      verbose_( parameters_->verbose() ),
      initialized_( false )
    {}

    double time () const { return timeProvider_.time(); }
    double timeStepSize () const { return timeProvider_.deltaT(); }

    void initialTimeStepSize ( double explicitEstimate )
    {
      timeProvider_.provideTimeStepEstimate( parameters().initialDeltaT( explicitEstimate ) );
      initialized_ = true;
    }

    template< class Monitor >
    void reduceTimeStep ( double explicitEstimate, const Monitor &monitor )
    {
      if( !initialized_ )
        DUNE_THROW( Dune::InvalidStateException, "ImplicitRungeKuttaSolver must be initialized before first solve." );

      double factor( 1 );
      parameters().cflFactor( explicitEstimate, explicitEstimate, monitor.newtonIterations_, false, factor );

      if( !((factor >= std::numeric_limits< double >::min()) && (factor < 1.0)) )
        DUNE_THROW( Dune::InvalidStateException, "invalid cfl factor: " << factor );

      cfl_ *= factor;
      timeProvider_.provideTimeStepEstimate( cfl_ * explicitEstimate );

      if( (verbose_ >= ImplicitRungeKuttaSolverParameters::cflVerbosity) && (Dune::Fem::MPIManager::rank() == 0) )
      {
        Dune::derr << "Implicit Runge-Kutta step failed to convergence." << std::endl;
        Dune::derr << "New cfl number: " << cfl_
                   << ", iterations per time step: ILS = " << monitor.linearSolverIterations_
                   << ", INLS = " << monitor.newtonIterations_
                   << std::endl;
      }

      timeProvider_.provideTimeStepEstimate( factor * timeStepSize() );
      timeProvider_.invalidateTimeStep();
    }

    template< class Monitor >
    void timeStepEstimate ( double explicitEstimate, const Monitor &monitor )
    {
      if( !initialized_ )
        DUNE_THROW( Dune::InvalidStateException, "ImplicitRungeKuttaSolver must be initialized before first solve." );

      double factor( 1 );
      parameters().cflFactor( explicitEstimate, explicitEstimate, monitor.newtonIterations_, true, factor );
      if( !((factor >= std::numeric_limits< double >::min()) && (factor <= std::numeric_limits< double >::max())) )
        DUNE_THROW( Dune::InvalidStateException, "invalid cfl factor: " << factor );

      const double oldCfl = cfl_;
      cfl_ = std::min( cflMax_, factor * cfl_ );
      timeProvider_.provideTimeStepEstimate( cfl_ * explicitEstimate );

      if( (cfl_ != oldCfl) && (verbose_ >= ImplicitRungeKuttaSolverParameters::cflVerbosity) && (Dune::Fem::MPIManager::rank() == 0) )
      {
        Dune::derr << "New cfl number: " << cfl_
                   << ", iterations per time step: ILS = " << monitor.linearSolverIterations_
                   << ", INLS = " << monitor.newtonIterations_
                   << std::endl;
      }
    }

  private:
    const ParametersType &parameters () const
    {
      assert( parameters_ );
      return *parameters_;
    }

    TimeProviderType &timeProvider_;
    const ParametersType *parameters_;
    double cfl_, cflMax_;
    int verbose_;
    bool initialized_;
  };



  /** \brief Implicit RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver >
  class ImplicitRungeKuttaSolver
  : public BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, ImplicitRungeKuttaTimeStepControl >
  {
    typedef ImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver > ThisType;
    typedef BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, ImplicitRungeKuttaTimeStepControl > BaseType;

  public:
    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef typename BaseType::TimeStepControlType TimeStepControlType;


    typedef typename TimeStepControlType::TimeProviderType TimeProviderType;
    typedef typename TimeStepControlType::ParametersType ParametersType;

    /** \brief constructor 
     *
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  butcherTable  butcher table to use
     *  \param[in]  traits        additional traits
     */
    ImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                               TimeProviderType &timeProvider,
                               int order = 1,
                               const ParametersType &parameters = ParametersType() )
    : BaseType( helmholtzOp, butcherTable( order ), TimeStepControlType( timeProvider, parameters ) )
    {}

  protected:
    static SimpleButcherTable< double > butcherTable ( int order )
    {
      switch( order )
      {
      case 1:
        return implicitEulerButcherTable();
      case 2:
        return gauss2ButcherTable();
      case 3:
        return implicit3ButcherTable();
      case 4:
        return implicit34ButcherTable();
      default:
        DUNE_THROW( NotImplemented, "Implicit Runge-Kutta method of order " << order << " not implemented." );
      }
    }
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH
