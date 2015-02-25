#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_TIMESTEPCONTROL_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_TIMESTEPCONTROL_HH

//- system includes
#include <cassert>
#include <memory>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/io/parameter.hh>
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

  protected:
    // number of minimal iterations that the linear solver should do
    // if the number of iterations done is smaller then the cfl number is increased
    const int minIter_;
    // number of maximal iterations that the linear solver should do
    // if the number of iterations larger then the cfl number is decreased
    const int maxIter_;
    // factor for cfl number on increase (decrease is 0.5)
    const double sigma_;

  public:
    ImplicitRungeKuttaSolverParameters ()
    : minIter_( Dune::Fem::Parameter::getValue< int >( "fem.ode.miniterations" , 14 ) ),
      maxIter_( Dune::Fem::Parameter::getValue< int >( "fem.ode.maxiterations" , 16 ) ),
      sigma_( Dune::Fem::Parameter::getValue< double >( "fem.ode.cflincrease" , 1.1 ) )
    {}

    // destructor (virtual)
    virtual ~ImplicitRungeKuttaSolverParameters() {}

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
        if (iter < minIter_)
        {
          if( imOpTimeStepEstimate <= exOpTimeStepEstimate )
          {
            factor = sigma_;
            changed = true;
          }
        }
        else if (iter > maxIter_)
        {
          factor = (double)maxIter_/(sigma_*(double)iter);
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

    // return number of max linear iterations per newton step
    virtual int maxLinearIterations () const { return maxIter_; }
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

    void initialTimeStepSize ( double helmholtzEstimate, double sourceTermEstimate )
    {
      double estimate = std::numeric_limits< double >::max();
      parameters().initTimeStepEstimate( sourceTermEstimate, helmholtzEstimate, estimate, cfl_ );
      timeProvider_.provideTimeStepEstimate( estimate );
      initialized_ = true;
    }

    template< class Monitor >
    void reduceTimeStep ( double helmholtzEstimate, double sourceTermEstimate, const Monitor &monitor )
    {
      if( !initialized_ )
        DUNE_THROW( Dune::InvalidStateException, "ImplicitRungeKuttaSolver must be initialized before first solve." );

      double factor( 1 );
      parameters().cflFactor( helmholtzEstimate, sourceTermEstimate, monitor.linearSolverIterations_, false, factor );

      if( !((factor >= std::numeric_limits< double >::min()) && (factor < 1.0)) )
        DUNE_THROW( Dune::InvalidStateException, "invalid cfl factor: " << factor );

      cfl_ *= factor;

      if( (verbose_ >= ImplicitRungeKuttaSolverParameters::noConvergenceVerbosity) && (Dune::Fem::MPIManager::rank() == 0) )
      {
        Dune::derr << "Implicit Runge-Kutta step failed to converge." << std::endl;
        Dune::derr << "New cfl number: " << cfl_
                   << ", iterations per time step: ILS = " << monitor.linearSolverIterations_
                   << ", INLS = " << monitor.newtonIterations_
                   << std::endl;
      }

      timeProvider_.provideTimeStepEstimate( factor * timeStepSize() );
      timeProvider_.invalidateTimeStep();
    }

    template< class Monitor >
    void timeStepEstimate ( double helmholtzEstimate, double sourceTermEstimate, const Monitor &monitor )
    {
      if( !initialized_ )
        DUNE_THROW( Dune::InvalidStateException, "ImplicitRungeKuttaSolver must be initialized before first solve." );

      double factor( 1 );
      // true means converged, which is always true since this function is only called
      // when the implicit solver did converge
      parameters().cflFactor( helmholtzEstimate, sourceTermEstimate, monitor.linearSolverIterations_, true, factor );
      if( !((factor >= std::numeric_limits< double >::min()) && (factor <= std::numeric_limits< double >::max())) )
        DUNE_THROW( Dune::InvalidStateException, "invalid cfl factor: " << factor );

      const double oldCfl = cfl_;
      cfl_ = std::min( cflMax_, factor * cfl_ );

      timeProvider_.provideTimeStepEstimate( std::min( sourceTermEstimate, cfl_ * helmholtzEstimate ) );

      if( (cfl_ != oldCfl) && (verbose_ >= ImplicitRungeKuttaSolverParameters::cflVerbosity) && (Dune::Fem::MPIManager::rank() == 0) )
      {
        Dune::derr << "New cfl number: " << cfl_
                   << ", iterations per time step: ILS = " << monitor.linearSolverIterations_
                   << ", INLS = " << monitor.newtonIterations_
                   << std::endl;
      }
    }

    bool computeError () const { return false; }

  protected:
    const ParametersType &parameters () const
    {
      assert( parameters_ );
      return *parameters_;
    }

    TimeProviderType &timeProvider_;
    std::shared_ptr< const ParametersType > parameters_;
    double cfl_, cflMax_;
    int verbose_;
    bool initialized_;
  };



  // PIDTimeStepControl
  // ------------------

  /** \brief PID time step control

      See also:
        D. Kuzmin and S.Turek. Numerical simulation of turbulent bubbly flows. Techreport Uni Dortmund. 2004

      and the original article:
        Valli, Coutinho, and Carey. Adaptive Control for Time Step Selection in Finite Element
        Simulation of Coupled Viscous Flow and Heat Transfer. Proc of the 10th
        International Conference on Numerical Methods in Fluids. 1998.
   */
  class PIDTimeStepControl : public ImplicitRungeKuttaTimeStepControl
  {
    typedef PIDTimeStepControl ThisType;
    typedef ImplicitRungeKuttaTimeStepControl BaseType;

  protected:
    using BaseType :: initialized_;
    using BaseType :: cfl_;
    using BaseType :: parameters ;
  public:
    typedef Dune::Fem::TimeProviderBase TimeProviderType;
    typedef ImplicitRungeKuttaSolverParameters ParametersType;

    explicit PIDTimeStepControl ( TimeProviderType &timeProvider,
                                  const ParametersType &parameters = ParametersType() )
    : BaseType( timeProvider, parameters ),
      errors_(),
      tol_( 1e-3 )
    {
      if( Dune::Fem::Parameter::getValue("fem.ode.pidcontrol", bool(false) ) )
      {
        tol_ = Dune::Fem::Parameter::getValue("fem.ode.pidtolerance", tol_ );
        errors_.resize( 3, tol_ );
      }
    }

    bool computeError () const { return ! errors_.empty() ; }

    template< class Monitor >
    void timeStepEstimate ( double helmholtzEstimate, double sourceTermEstimate, const Monitor &monitor )
    {
      if( !initialized_ )
        DUNE_THROW( Dune::InvalidStateException, "ImplicitRungeKuttaSolver must be initialized before first solve." );

      if( computeError() ) // use pid control
      {
        cfl_ = 1.0; // reset cfl for next reduceTimeStep
        double dtEst = pidTimeStepControl( std::min( sourceTermEstimate, helmholtzEstimate ), monitor );
        const int targetIterations = parameters().maxLinearIterations();
        /*
        if( monitor.linearSolverIterations_ > targetIterations &&
            targetIterations > 0 )
        {
          dtEst *= double( targetIterations ) / double(monitor.linearSolverIterations_);
        }
        */
        std::cout << "Set dt = " << dtEst << std::endl;
        timeProvider_.provideTimeStepEstimate( dtEst );

        if( (verbose_ >= ImplicitRungeKuttaSolverParameters::cflVerbosity) && (Dune::Fem::MPIManager::rank() == 0) )
        {
          Dune::derr << "New dt: " << dtEst
                     << ", iterations per time step: ILS = " << monitor.linearSolverIterations_
                     << ", INLS = " << monitor.newtonIterations_
                     << std::endl;
        }
      }
      else
      {
        BaseType::timeStepEstimate( helmholtzEstimate, sourceTermEstimate, monitor );
      }
    }

    template < class Monitor >
    double pidTimeStepControl( const double dt, const Monitor& monitor )
    {
      // get error || u^n - u^n+1 || / || u^n+1 || from monitor
      const double error = monitor.error_;
      std::cout << error << " error " << std::endl;
      if( std::abs( error ) < 1e-12 ) return 10. * dt;

      // shift errors
      for( int i=0; i<2; ++i )
      {
        errors_[ i ] = errors_[i+1];
      }

      // store new error
      errors_[ 2 ] = error ;

      if( error > tol_ )
      {
        // adjust dt by given tolerance
        const double newDt = dt * tol_ / error;
        return newDt;
      }
      else if( error > 1e-12 )
      {
        // values taking from turek time stepping paper
        const double kP = 0.075 ;
        const double kI = 0.175 ;
        const double kD = 0.01 ;
        const double newDt = (dt * std::pow( errors_[ 1 ] / errors_[ 2 ], kP ) *
                             std::pow( tol_         / errors_[ 2 ], kI ) *
                             std::pow( errors_[0]*errors_[0]/errors_[ 1 ]/errors_[ 2 ], kD ));
        std::cout << "newDt = " << newDt << std::endl;
        return newDt;
      }

      return dt ;
    }

  protected:
    std::vector< double > errors_;
    double tol_;
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_TIMESTEPCONTROL_HH
