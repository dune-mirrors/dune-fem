#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICROW_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_BASICROW_HH

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

  // NoROWRungeKuttaSourceTerm
  // ------------------------------

  struct NoROWRungeKuttaSourceTerm
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

  struct ROWSolverParameter
#ifndef DOXYGEN
  : public LocalParameter< ROWSolverParameter, ROWSolverParameter >
#endif
  {
    ROWSolverParameter () {}

    virtual double linAbsTolParameter ( )  const
    {
      return Parameter::getValue< double >( "fem.solver.row.linabstol", 1e-6 );
    }

    virtual double linReductionParameter ( ) const
    {
      return Parameter::getValue< double >( "fem.solver.row.linreduction", 1e-4  );
    }

    virtual bool linearSolverVerbose () const
    {
      return Parameter::getValue< bool >( "fem.solver.row.linear.verbose", false );
    }

    virtual int maxLinearIterationsParameter () const
    {
      return Parameter::getValue< int >( "fem.solver.row.maxlineariterations", std::numeric_limits< int >::max() );
    }
  };



  /** \brief ROW RungeKutta ODE solver. */
  template< class HelmholtzOperator, class NonlinearSolver, class TimeStepControl, class SourceTerm = NoROWRungeKuttaSourceTerm >
  class BasicROWRungeKuttaSolver
  : public OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType >
  {
    typedef BasicROWRungeKuttaSolver< HelmholtzOperator, NonlinearSolver, TimeStepControl, SourceTerm > ThisType;
    typedef OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType > BaseType;

  public:
    typedef typename BaseType::MonitorType MonitorType;
    typedef typename BaseType::DestinationType DestinationType;

    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef NonlinearSolver NonlinearSolverType;
    typedef TimeStepControl TimeStepControlType;
    typedef SourceTerm SourceTermType;

    typedef typename HelmholtzOperator::SpaceOperatorType::PreconditionOperatorType PreconditionOperatorType;

    typedef Dune::Fem::TimeProviderBase TimeProviderType;

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp      Helmholtz operator \f$L\f$
     *  \param[in]  butcherTable     butcher table to use
     *  \param[in]  timeStepControl  time step controller
     *  \param[in]  sourceTerm       additional source term
     */
    template< class ButcherTable >
    BasicROWRungeKuttaSolver ( HelmholtzOperatorType &helmholtzOp,
                                    const ButcherTable &butcherTable,
                                    const TimeStepControlType &timeStepControl = TimeStepControl(),
                                    const SourceTermType &sourceTerm = SourceTermType(),
                                    const ROWSolverParameter &parameter = ROWSolverParameter() )
    : helmholtzOp_( helmholtzOp ),
      nonlinearSolver_( helmholtzOp_ ),
      timeStepControl_( timeStepControl ),
      sourceTerm_( sourceTerm ),
      stages_( butcherTable.stages() ),
      alpha_( butcherTable.A() ),
      alpha2_( butcherTable.B() ),
      gamma_( stages() ),
      beta_( stages() ),
      c_( butcherTable.c() ),
      rhs_( "RK rhs", helmholtzOp_.space() ),
      temp_( "RK temp", helmholtzOp_.space() ),
      update_( stages(), nullptr ),
      linAbsTol_( parameter.linAbsTolParameter( ) ),
      linReduction_( parameter.linReductionParameter( ) ),
      linVerbose_( parameter.linearSolverVerbose() ),
      maxLinearIterations_( parameter.maxLinearIterationsParameter() ),
      preconditioner_(helmholtzOp.spaceOperator().preconditioner())
    {
      std::cout << "ROW method of order=" << butcherTable.order() << " with " << stages_ << " stages" << std::endl;
      // create intermediate functions
      for( int i = 0; i < stages(); ++i )
      {
        std::ostringstream name;
        name << "RK stage " << i;
        update_[ i ] = new DestinationType( name.str(), helmholtzOp_.space() );
      }

      // compute coefficients
      for( int i = 0; i < stages(); ++i )
        gamma_[ i ] = alpha_[ i ][ i ];

      alpha_.invert();
      alpha_.mtv( butcherTable.b(), beta_ );
      alpha2_.rightmultiply( alpha_ );
    }

    /** \brief destructor */
    ~BasicROWRungeKuttaSolver ()
    {
      for( int i = 0; i < stages(); ++i )
        delete update_[ i ];
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

      typename HelmholtzOperatorType::JacobianOperatorType jOp( "jacobianOperator", U.space(), U.space() );

      const double time = timeStepControl_.time();
      const double timeStepSize = timeStepControl_.timeStepSize();
      assert( timeStepSize > 0.0 );
      // std::cout << "solving... " << time << "    :     " << U << std::endl;
      for( int s = 0; s < stages(); ++s )
      {
        // update for stage s
        DestinationType& updateStage = *update_[ s ];

        rhs_.assign( U );
        for( int k = 0; k < s; ++k )
          rhs_.axpy( alpha2_[ s ][ k ], *update_[ k ] );
        helmholtzOp_.spaceOperator()(rhs_,updateStage);
        updateStage *= (gamma_[s]*timeStepSize);
        for( int k = 0; k < s; ++k )
          updateStage.axpy( -gamma_[s]*alpha_[ s ][ k ], *update_[ k ] );

        rhs_.assign( updateStage );

        // solve the system
        const double stageTime = time + c_[ s ]*timeStepSize;
        helmholtzOp_.setTime( stageTime );
        // compute jacobian if the diagonal entry in the butcher tableau changes
        // if ( s==0 || (gamma_[s-1] != gamma_[s]) )
        {
          // std::cout << "lambda=" << gamma_[ s ]*timeStepSize << std::endl;
          helmholtzOp_.setLambda( gamma_[ s ]*timeStepSize );
          helmholtzOp_.jacobian( U, jOp );
        }
        const int remLinearIts = maxLinearIterations_;
        if (preconditioner_)
        {
          typename NonlinearSolverType::LinearInverseOperatorType jInv( jOp, *preconditioner_, linReduction_, linAbsTol_, remLinearIts, linVerbose_ );
          jInv( rhs_, updateStage );
          monitor.linearSolverIterations_ += jInv.iterations();
        }
        else
        {
          typename NonlinearSolverType::LinearInverseOperatorType jInv( jOp, linReduction_, linAbsTol_, remLinearIts, linVerbose_ );
          jInv( rhs_, updateStage );
          monitor.linearSolverIterations_ += jInv.iterations();
        }
      }

      double error = 0.0;
      if(0 && timeStepControl_.computeError() )
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

    HelmholtzOperatorType &helmholtzOp_;
    NonlinearSolverType nonlinearSolver_;
    TimeStepControl timeStepControl_;
    SourceTerm sourceTerm_;

    int stages_;
    double delta_;
    Dune::DynamicMatrix< double > alpha_, alpha2_;
    Dune::DynamicVector< double > gamma_, beta_, c_;

    DestinationType rhs_,temp_;
    std::vector< DestinationType * > update_;

    const double linAbsTol_, linReduction_;
    const bool linVerbose_;
    const int maxLinearIterations_;

    const PreconditionOperatorType *preconditioner_;
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICROW_HH
