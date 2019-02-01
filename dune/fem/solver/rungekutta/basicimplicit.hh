#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH

//- system includes
#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

//- dune-fem includes
#include <dune/fem/solver/odesolverinterface.hh>

namespace DuneODE
{

  // NoImplicitRungeKuttaSourceTerm
  // ------------------------------

  struct NoImplicitRungeKuttaSourceTerm
  {
    template< class T >
    bool operator() ( double time, double timeStepSize, int stage, const T &u, const std::vector< T * > &update, T &source )
    {
      return false;
    }

    template< class T >
    void limit( T& update, const double time ) {}

    template< class T >
    double initialTimeStepEstimate ( double time, const T &u ) const
    {
      // return negative value to indicate that implicit time step should be used
      return -1.0;
    }

    double timeStepEstimate () const
    {
      return std::numeric_limits< double >::max();
    }
  };



  /** \brief Implicit RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver, class TimeStepControl, class SourceTerm = NoImplicitRungeKuttaSourceTerm >
  class BasicImplicitRungeKuttaSolver
  : public OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType >
  {
    typedef BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl, SourceTerm > ThisType;
    typedef OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType > BaseType;

  public:
    typedef typename BaseType::MonitorType MonitorType;
    typedef typename BaseType::DestinationType DestinationType;

    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef NonlinearSolver NonlinearSolverType;
    typedef TimeStepControl TimeStepControlType;
    typedef SourceTerm SourceTermType;

    typedef Dune::Fem::TimeProviderBase TimeProviderType;

    typedef typename TimeStepControlType::ParameterType ParameterType;
    typedef typename NonlinearSolver::ParameterType     NonlinearSolverParameterType;

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp      Helmholtz operator \f$L\f$
     *  \param[in]  butcherTable     butcher table to use
     *  \param[in]  timeStepControl  time step controller
     *  \param[in]  sourceTerm       additional source term
     */
    template< class ButcherTable >
    BasicImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                                    const ButcherTable &butcherTable,
                                    const TimeStepControlType &timeStepControl,
                                    const SourceTermType &sourceTerm,
                                    const NonlinearSolverParameterType& parameters )
    : helmholtzOp_( helmholtzOp ),
      nonlinearSolver_( parameters ),
      timeStepControl_( timeStepControl ),
      sourceTerm_( sourceTerm ),
      stages_( butcherTable.stages() ),
      alpha_( butcherTable.A() ),
      gamma_( stages() ),
      beta_( stages() ),
      c_( butcherTable.c() ),
      rhs_( "RK rhs", helmholtzOp_.space() ),
      updateStorage_(),
      update_( stages(), nullptr )
    {
      setup( butcherTable );
    }

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp      Helmholtz operator \f$L\f$
     *  \param[in]  butcherTable     butcher table to use
     *  \param[in]  timeStepControl  time step controller
     */
    template< class ButcherTable >
    BasicImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                                    const ButcherTable &butcherTable,
                                    const TimeStepControlType &timeStepControl,
                                    const NonlinearSolverParameterType& parameters )
    : helmholtzOp_( helmholtzOp ),
      nonlinearSolver_( parameters ),
      timeStepControl_( timeStepControl ),
      stages_( butcherTable.stages() ),
      alpha_( butcherTable.A() ),
      gamma_( stages() ),
      beta_( stages() ),
      c_( butcherTable.c() ),
      rhs_( "RK rhs", helmholtzOp_.space() ),
      updateStorage_(),
      update_( stages(), nullptr )
    {
      setup( butcherTable );
    }

    template< class ButcherTable >
    BasicImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                                    const ButcherTable &butcherTable,
                                    const TimeStepControlType &timeStepControl,
                                    const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BasicImplicitRungeKuttaSolver( helmholtzOp, butcherTable, timeStepControl,
        NonlinearSolverParameterType( parameter ) )
    {
    }

    template< class ButcherTable >
    BasicImplicitRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                                    const ButcherTable &butcherTable,
                                    const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BasicImplicitRungeKuttaSolver( helmholtzOp, butcherTable,
        TimeStepControlType( parameter ), NonlinearSolverParameterType( parameter ) )
    {
    }

    template< class ButcherTable >
    void setup( const ButcherTable& butcherTable )
    {
      update_.clear();
      update_.resize( stages(), nullptr );
      updateStorage_.resize( stages() );

      // create intermediate functions
      for( int i = 0; i < stages(); ++i )
      {
        std::ostringstream name;
        name << "RK stage " << i;
        updateStorage_[ i ].reset( new DestinationType( name.str(), helmholtzOp_.space() ) );
        update_[ i ] = updateStorage_[ i ].operator ->();
      }

      // compute coefficients
      Dune::DynamicMatrix< double > AL( alpha_ );
      for( int i = 0; i < stages(); ++i )
      {
        gamma_[ i ] = AL[ i ][ i ];
        AL[ i ][ i ] = 0.0;
      }

      alpha_.invert();
      alpha_.mtv( butcherTable.b(), beta_ );

      alpha_.leftmultiply( AL );
      for( int i = 0; i < stages(); ++i )
        alpha_[ i ][ i ] = gamma_[ i ];

      for( int i = 0; i < stages(); ++i )
      {
        gamma_[ i ] = 1.0;
        for( int j = 0; j < i; ++j )
          gamma_[ i ] -= alpha_[ i ][ j ];
      }

      delta_ = 1.0;
      for( int i = 0; i < stages(); ++i )
        delta_ -= beta_[ i ];

    }

    //! apply operator once to get dt estimate
    void initialize ( const DestinationType &U0 )
    {
      const double time = timeStepControl_.time();

      helmholtzOp_.setTime( time );
      helmholtzOp_.initializeTimeStepSize( U0 );
      const double helmholtzEstimate = helmholtzOp_.timeStepEstimate();

      double sourceTermEstimate = sourceTerm_.initialTimeStepEstimate( time, U0 );
      // negative time step is given by the empty source term
      if( sourceTermEstimate < 0.0 ) sourceTermEstimate = helmholtzEstimate ;

      timeStepControl_.initialTimeStepSize( helmholtzEstimate, sourceTermEstimate );
    }

    using BaseType::solve;

    //! solve the system
    void solve ( DestinationType &U, MonitorType &monitor )
    {
      monitor.reset();

      const double time = timeStepControl_.time();
      const double timeStepSize = timeStepControl_.timeStepSize();
      assert( timeStepSize > 0.0 );
      for( int s = 0; s < stages(); ++s )
      {
        assert( update_[ s ] );
        // update for stage s
        DestinationType& updateStage = *update_[ s ];

        // assemble rhs of nonlinear equation
        updateStage.assign( U );
        updateStage *= gamma_[ s ];
        for( int k = 0; k < s; ++k )
          updateStage.axpy( alpha_[ s ][ k ], *update_[ k ] );

        const double stageTime = time + c_[ s ]*timeStepSize;
        if( sourceTerm_( time, timeStepSize, s, U, update_, rhs_ ) )
        {
          updateStage.axpy( alpha_[ s ][ s ]*timeStepSize, rhs_ );
          sourceTerm_.limit( updateStage, stageTime );
        }

        // apply Helmholtz operator to right hand side
        helmholtzOp_.setTime( stageTime );
        helmholtzOp_.setLambda( 0 );
        helmholtzOp_( updateStage, rhs_ );

        // solve the system
        helmholtzOp_.setLambda( alpha_[ s ][ s ]*timeStepSize );
        nonlinearSolver_.bind( helmholtzOp_ );
        nonlinearSolver_( rhs_, updateStage );
        nonlinearSolver_.unbind();

        // update monitor
        monitor.newtonIterations_       += nonlinearSolver_.iterations();
        monitor.linearSolverIterations_ += nonlinearSolver_.linearIterations();

        // on failure break solving
        if( !nonlinearSolver_.converged() )
          return timeStepControl_.reduceTimeStep( helmholtzOp_.timeStepEstimate(), sourceTerm_.timeStepEstimate(), monitor );
      }

      double error = 0.0;
      if( timeStepControl_.computeError() )
      {
        // store U (to be revised)
        DestinationType Uerr( U );

        // update solution
        U *= delta_;
        for( int s = 0; s < stages(); ++s )
          U.axpy( beta_[ s ], *update_[ s ] );

        //error = infNorm( U, Uerr );
        Uerr.axpy( -1.0, U );
        const double errorU = Uerr.scalarProductDofs( Uerr );
        const double normU = U.scalarProductDofs( U );

        if( normU > 0 && errorU > 0 )
        {
          error = std::sqrt( errorU / normU );
        }
        std::cout << std::scientific << "Error in RK = " << error << " norm " << errorU << " " << normU << std::endl;
        //std::cout << std::scientific << "Error in RK = " << error << std::endl;
      }
      else
      {
        // update solution
        U *= delta_;
        for( int s = 0; s < stages(); ++s )
          U.axpy( beta_[ s ], *update_[ s ] );
      }
      // set error to monitor
      monitor.error_ = error;

      // update time step size
      timeStepControl_.timeStepEstimate( helmholtzOp_.timeStepEstimate(), sourceTerm_.timeStepEstimate(), monitor );
    }

    int stages () const { return stages_; }

    void description ( std::ostream &out ) const
    {
      out << "Generic " << stages() << "-stage implicit Runge-Kutta solver.\\\\" << std::endl;
    }

  protected:
    double infNorm(const DestinationType& U, const DestinationType& Uerr ) const
    {
      typedef typename DestinationType :: ConstDofIteratorType ConstDofIteratorType ;
      const ConstDofIteratorType uend = U.dend();
      double res = 0;
      for( ConstDofIteratorType u = U.dbegin(), uerr = Uerr.dbegin(); u != uend; ++u, ++uerr )
      {
        double uval = *u;
        double uerrval = *uerr ;
        double div = std::abs( std::max( uval, uerrval ) );

        double norm = std::abs( uval - uerrval );
        if( std::abs(div) > 1e-12 )
          norm /= div;
        res = std::max( res, norm );
      }
      return res;
    }

    HelmholtzOperatorType&   helmholtzOp_;
    NonlinearSolverType      nonlinearSolver_;
    TimeStepControl          timeStepControl_;
    SourceTerm               sourceTerm_;

    int stages_;
    double delta_;
    Dune::DynamicMatrix< double > alpha_;
    Dune::DynamicVector< double > gamma_, beta_, c_;

    DestinationType rhs_;
    std::vector< std::unique_ptr< DestinationType > > updateStorage_;
    std::vector< DestinationType* > update_;
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH
