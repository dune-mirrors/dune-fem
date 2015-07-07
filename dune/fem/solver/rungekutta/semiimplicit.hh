#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_SEMIIMPLICIT_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_SEMIIMPLICIT_HH

//- system includes
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-fem includes
#include <dune/fem/solver/rungekutta/basicimplicit.hh>
#include <dune/fem/solver/rungekutta/butchertable.hh>
#include <dune/fem/solver/rungekutta/timestepcontrol.hh>

namespace DuneODE
{

  // SemiImplicitRungeKuttaSourceTerm
  // --------------------------------

  template< class ExplicitOperator >
  class SemiImplicitRungeKuttaSourceTerm
  {
    typedef SemiImplicitRungeKuttaSourceTerm< ExplicitOperator > ThisType;

  public:
    typedef ExplicitOperator ExplicitOperatorType;

    typedef typename ExplicitOperatorType::DestinationType DestinationType;

    template< class ButcherTable >
    SemiImplicitRungeKuttaSourceTerm ( ExplicitOperatorType &explicitOp,
                                       const ButcherTable &butcherTable,
                                       const Dune::DynamicMatrix< double > &implicitA )
    : explicitOp_( explicitOp ),
      alpha_( butcherTable.A() ),
      gamma_( butcherTable.stages() ),
      c_( butcherTable.c() ),
      uex_( "SIRK u-explicit", explicitOp_.space() ),
      limiter_( explicitOp_.hasLimiter() )
    {
      Dune::DynamicMatrix< double > Ainv( implicitA );
      Ainv.invert();
      alpha_.rightmultiply( Ainv );

      for( int i = 0; i < butcherTable.stages(); ++i )
      {
        gamma_[ i ] = 1.0;
        for( int j = 0; j < i; ++j )
          gamma_[ i ] -= alpha_[ i ][ j ];
      }
    }

    bool operator() ( double time, double timeStepSize, int stage,
                      const DestinationType &u, const std::vector< DestinationType * > &update,
                      DestinationType &source )
    {
      uex_.assign( u );
      uex_ *= gamma_[ stage ];
      for( int k = 0; k < stage; ++k )
        uex_.axpy( alpha_[ stage ][ k ], *update[ k ] );
      explicitOp_.setTime( time + c_[ stage ]*timeStepSize );
      explicitOp_( uex_, source );

      return true;
    }

    void limit( DestinationType& update, const double time )
    {
      if( limiter_ )
      {
        // set correct time
        explicitOp_.setTime( time );
        // copy given function
        uex_.assign( update );
        // apply limiter
        explicitOp_.limit( uex_, update );
      }
    }

    double initialTimeStepEstimate ( double time, const DestinationType &u ) const
    {
      explicitOp_.setTime( time );
      explicitOp_.initializeTimeStepSize( u );
      return explicitOp_.timeStepEstimate();
    }

    double timeStepEstimate () const
    {
      return explicitOp_.timeStepEstimate();
    }

  private:
    ExplicitOperatorType &explicitOp_;
    Dune::DynamicMatrix< double > alpha_;
    Dune::DynamicVector< double > gamma_, c_;
    DestinationType uex_;
    const bool limiter_ ;
  };



  // SemiImplicitRungeKuttaSolver
  // ----------------------------

  /** \brief Implicit RungeKutta ODE solver. */
  template< class ExplicitOperator, class HelmholtzOperator, class NonlinearSolver >
  class SemiImplicitRungeKuttaSolver
  : public BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, ImplicitRungeKuttaTimeStepControl, SemiImplicitRungeKuttaSourceTerm< ExplicitOperator > >
  {
    typedef SemiImplicitRungeKuttaSolver< ExplicitOperator, HelmholtzOperator, NonlinearSolver > ThisType;
    typedef BasicImplicitRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, ImplicitRungeKuttaTimeStepControl, SemiImplicitRungeKuttaSourceTerm< ExplicitOperator > > BaseType;

  public:
    typedef ExplicitOperator ExplicitOperatorType;
    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef typename BaseType::TimeStepControlType TimeStepControlType;
    typedef typename BaseType::SourceTermType SourceTermType;

    typedef typename TimeStepControlType::TimeProviderType TimeProviderType;
    typedef typename BaseType::ParametersType                  ParametersType;
    typedef typename BaseType::NonlinearSolverParametersType   NonlinearSolverParametersType;

    /** \brief constructor
     *
     *  \param[in]  explicitOp    explicit operator
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  order         order of butcher table to use
     *  \param[in]  tscParam      parameters for implicit time step control
     *  \param[in]  nlsParam      parameters for non linear solver control
     */
    SemiImplicitRungeKuttaSolver ( ExplicitOperatorType &explicitOp,
                                   HelmholtzOperatorType &helmholtzOp,
                                   TimeProviderType &timeProvider,
                                   int order = 1,
                                   const ParametersType& tscParams = ParametersType(),
                                   const NonlinearSolverParametersType& nlsParams = NonlinearSolverParametersType() )
    : BaseType( helmholtzOp,
                butcherTable( order, false ),
                TimeStepControlType( timeProvider, tscParams ),
                SourceTermType( explicitOp, butcherTable( order, true ), butcherTable( order, false ).A() ),
                nlsParams )
    {}

    /** \brief constructor
     *
     *  \param[in]  explicitOp    explicit operator
     *  \param[in]  helmholtzOp   Helmholtz operator \f$L\f$
     *  \param[in]  timeProvider  time provider
     *  \param[in]  tscParam      parameters for implicit time step control
     *  \param[in]  nlsParam      parameters for non linear solver control
     */
    SemiImplicitRungeKuttaSolver ( ExplicitOperatorType &explicitOp,
                                   HelmholtzOperatorType &helmholtzOp,
                                   TimeProviderType &timeProvider,
                                   const ParametersType& tscParams = ParametersType(),
                                   const NonlinearSolverParametersType& nlsParams = NonlinearSolverParametersType() )
    : BaseType( helmholtzOp,
                butcherTable( -1, false ),
                TimeStepControlType( timeProvider, tscParams ),
                SourceTermType( explicitOp, butcherTable( -1, true ), butcherTable( -1, false ).A() ),
                nlsParams )
    {}
  protected:
    static SimpleButcherTable< double > butcherTable ( int order, bool expl )
    {
      switch( order )
      {
      case 1:
        return semiImplicitEulerButcherTable( expl );
      case 2:
        return semiImplicitSSP222ButcherTable( expl );
      case 3:
        return semiImplicit33ButcherTable( expl );
      case 4:
        return semiImplicitIERK45ButcherTable( expl );
      default:
        DUNE_THROW( NotImplemented, "Semi-implicit Runge-Kutta method of order " << order << " not implemented." );
      }
    }
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_IMPLICIT_HH
