#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH

//- system includes 
#include <cassert>
#include <limits>
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/nullptr.hh>

//- dune-fem includes
#include <dune/fem/solver/odesolver.hh>

namespace DuneODE 
{

  // NoImplicitRungeKuttaSourceTerm
  // ------------------------------

  struct NoImplicitRungeKuttaSourceTerm
  {
    template< class T >
    bool source ( const T &u, const std::vector< T * > &update, int stage, T &source )
    {
      return false;
    }

    double timeStepEstimate ()
    {
      return std::numeric_limits< double >::max();
    }
  };



  /** \brief Implicit RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver, class TimeStepControl, class SourceTerm = NoImplicitRungeKuttaSourceTerm >
  class BasicImplicitRungeKuttaSolver
  : public OdeSolverInterface< typename HelmholtzOperator::DestinationType >
  {
    typedef BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl, SourceTerm > ThisType;
    typedef OdeSolverInterface< typename HelmholtzOperator::DestinationType > BaseType;

  public:
    typedef typename BaseType::MonitorType MonitorType;

    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef NonlinearSolver NonlinearSolverType;
    typedef TimeStepControl TimeStepControlType;
    typedef SourceTerm SourceTermType;

    typedef Dune::Fem::TimeProviderBase TimeProviderType;

    typedef typename HelmholtzOperatorType::DestinationType DestinationType;

    typedef typename DestinationType::DiscreteFunctionSpaceType SpaceType;
      
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
                                    const TimeStepControlType &timeStepControl = TimeStepControl(),
                                    const SourceTermType &sourceTerm = SourceTermType() )
    : helmholtzOp_( helmholtzOp ),
      nonlinearSolver_( helmholtzOp_ ),
      timeStepControl_( timeStepControl ),
      sourceTerm_( sourceTerm ),
      stages_( butcherTable.stages() ),
      alpha_( butcherTable.A() ),
      gamma_( stages() ),
      beta_( stages() ),
      c_( butcherTable.c() ),
      rhs_( "RK rhs", helmholtzOp_.space() ),
      update_( stages(), nullptr )
    {
      // create intermediate functions
      for( int i = 0; i < stages(); ++i )
      {
        std::ostringstream name;
        name << "RK stage " << i;
        update_[ i ] = new DestinationType( name.str(), helmholtzOp_.space() );
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

    /** \brief destructor */
    ~BasicImplicitRungeKuttaSolver ()
    {
      for( int i = 0; i < stages(); ++i )
        delete update_[ i ];
    }

    //! apply operator once to get dt estimate 
    void initialize ( const DestinationType &U0 )
    {
      helmholtzOp_( U0, *update_[ 0 ] );
      timeStepControl_.initialTimeStepSize( helmholtzOp_.timeStepEstimate() );
    }

    using BaseType::solve;

    //! solve the system 
    void solve ( DestinationType &U, MonitorType &monitor )
    {
      monitor.reset();

      const double time = timeStepControl_.time();
      const double timeStepSize = timeStepControl_.timeStepSize();
      assert( timeStepSize > 0.0 );

      for( int s; s < stages(); ++s )
      {
        // assemble rhs of nonlinear equation
        update_[ s ]->assign( U );
        *update_[ s ] *= gamma_[ s ];
        for( int k = 0; k < s; ++k )
          update_[ s ]->axpy( alpha_[ s ][ k ], *update_[ k ] );
        if( source( U, update_, s, rhs_ ) )
          update_[ s ]->axpy( alpha_[ s ][ s ]*timeStepSize, rhs_ ); 

        // apply Helmholtz operator to right hand side
        helmholtzOp_.setTime( time + c_[ s ]*timeStepSize );
        helmholtzOp_.setLambda( 0 );
        helmholtzOp_( *update_[ s ], rhs_ );

        // solve the system
        helmholtzOp_.setLambda( alpha_[ s ][ s ]*timeStepSize );
        nonlinearSolver_( rhs_, *update_[ s ] );

        // on failure break solving
        if( !nonlinearSolver_.converged() )
          return timeStepControl_.reduceTimeStep( helmholtzOp_.timeStepEstimate(), source.timeStepEstimate(), updateMonitor( monitor ) );
      }

      // update solution
      U *= delta_;
      for( int s = 0; s < stages(); ++s )
        U.axpy( beta_[ s ], *update_[ s ] );

      // update time step size
      timeStepControl_.timeStepEstimate( helmholtzOp_.timeStepEstimate(), source.timeStepEstimate(), updateMonitor( monitor ) );
    }

    int stages () const { return stages(); }

    void description ( std::ostream &out ) const
    {
      out << "Generic " << stages() << "-stage implicit Runge-Kutta solver.\\\\" << std::endl;
    }

  protected:
    const MonitorType &updateMonitor ( MonitorType &monitor ) const
    {
      monitor.newtonIterations_ = nonlinearSolver_.iterations();
      monitor.linearSolverIterations_ = nonlinearSolver_.linearIterations();
      return monitor;
    }

    HelmholtzOperatorType &helmholtzOp_;
    NonlinearSolverType nonlinearSolver_;
    TimeStepControl timeStepControl_;
    SourceTerm sourceTerm_;

    int stages_;
    double delta_;
    Dune::DynamicMatrix< double > alpha_;
    Dune::DynamicVector< double > gamma_, beta_, c_;

    DestinationType rhs_;
    std::vector< DestinationType * > update_;
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH
