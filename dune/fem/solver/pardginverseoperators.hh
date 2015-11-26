#ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/pardg.hh>

#ifdef USE_PARDG_ODE_SOLVER

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class DomainFunction, class RangeFunction = DomainFunction >
    class ParDGOperator;



    // ParDGOperator for AdaptiveDiscreteFunction
    // ------------------------------------------

    template< class DomainFunctionSpace, class RangeFunctionSpace >
    class ParDGOperator< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > >
    : public PARDG::Function
    {
      typedef ParDGOperator< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > > ThisType;

      typedef AdaptiveDiscreteFunction< DomainFunctionSpace > DomainFunctionType;
      typedef AdaptiveDiscreteFunction< RangeFunctionSpace > RangeFunctionType;

    public:
      typedef Operator< DomainFunctionType, RangeFunctionType > OperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      ParDGOperator ( const OperatorType &op, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
      : operator_( op ),
        domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace )
      {}

      void operator() ( const double *u, double *w, int i = 0 )
      {
        DomainFunctionType uFunction( "ParDGOperator u", domainSpace_, u );
        RangeFunctionType wFunction( "ParDGOperator w", rangeSpace_, w );
        operator_( uFunction, wFunction );
      }

      int dim_of_argument( int i = 0 ) const
      {
        assert( i == 0 );
        return domainSpace_.size();
      }

      int dim_of_value ( int i = 0 ) const
      {
        assert( i == 0 );
        return rangeSpace_.size();
      }

    private:
      const OperatorType &operator_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };




    // ParDGGeneralizedMinResInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    class ParDGGeneralizedMinResInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      typedef ParDGOperator< DiscreteFunction, DiscreteFunction > ParDGOperatorType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              double redEps, double absLimit,
                                              unsigned int maxIterations, bool verbose )
      : solver_( PARDG::Communicator::instance(), paramRestart() ),
        operator_( op ),
        preconditioner_( nullptr )
      {
        setupSolver( redEps, absLimit, maxIterations, verbose );
      }

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              double redEps, double absLimit,
                                              unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : solver_( PARDG::Communicator::instance(), paramRestart() ),
        operator_( op ),
        preconditioner_( nullptr )
      {
        setupSolver( redEps, absLimit, maxIterations, false );
      }

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              const PreconditionerType &preconditioner,
                                              double redEps, double absLimit,
                                              unsigned int maxIterations, bool verbose )
      : solver_( PARDG::Communicator::instance(), paramRestart() ),
        operator_( op ),
        preconditioner_( &preconditioner )
      {
        setupSolver( redEps, absLimit, maxIterations, verbose );
      }

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              const PreconditionerType &preconditioner,
                                              double redEps, double absLimit,
                                              unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : solver_( PARDG::Communicator::instance(), paramRestart() ),
        operator_( op ),
        preconditioner_( &preconditioner )
      {
        setupSolver( redEps, absLimit, maxIterations, false );
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ParDGOperatorType parDGOperator( operator_, w.space(), u.space() );
        if( preconditioner_ )
        {
          ParDGOperatorType parDGPreconditioner( *preconditioner_, w.space(), w.space() );
          solver_.set_preconditioner( parDGPreconditioner );
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
          solver_.unset_preconditioner();
        }
        else
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
      }

      unsigned int iterations () const
      {
        return solver_.number_of_iterations();
      }

    private:
      static int paramRestart ()
      {
        return Parameter::getValue< int >( "fem.solver.gmres.restart", 20 );
      }

      void setupSolver ( double redEps, double absLimit, unsigned int maxIterations, bool verbose )
      {
        PARDG::set_tolerance(solver_, redEps, absLimit, "fem.solver.errormeasure" );

        maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations );
        solver_.set_max_number_of_iterations( int( maxIterations ) );

        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose && Parameter :: verbose() )
        {
          solver_.IterativeSolver::set_output( std::cout );
          solver_.DynamicalObject::set_output( std::cout );
        }
      }

      mutable PARDG::GMRES solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
    };

    // ParDGBICGStabInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    class ParDGBiCGStabInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      typedef ParDGOperator< DiscreteFunction, DiscreteFunction > ParDGOperatorType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     double redEps, double absLimit,
                                     unsigned int maxIterations, bool verbose )
      : solver_( PARDG::Communicator::instance() ),
        operator_( op ),
        preconditioner_( nullptr )
      {
        setupSolver( redEps, absLimit, maxIterations, verbose );
      }

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     double redEps, double absLimit,
                                     unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : solver_( PARDG::Communicator::instance() ),
        operator_( op ),
        preconditioner_( nullptr )
      {
        setupSolver( redEps, absLimit, maxIterations, false );
      }

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     const PreconditionerType &preconditioner,
                                     double redEps, double absLimit,
                                     unsigned int maxIterations, bool verbose )
      : solver_( PARDG::Communicator::instance() ),
        operator_( op ),
        preconditioner_( &preconditioner )
      {
        setupSolver( redEps, absLimit, maxIterations, verbose );
      }

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     const PreconditionerType &preconditioner,
                                     double redEps, double absLimit,
                                     unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : solver_( PARDG::Communicator::instance() ),
        operator_( op ),
        preconditioner_( &preconditioner )
      {
        setupSolver( redEps, absLimit, maxIterations, false );
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ParDGOperatorType parDGOperator( operator_, w.space(), u.space() );
        if( preconditioner_ )
        {
          ParDGOperatorType parDGPreconditioner( *preconditioner_, w.space(), w.space() );
          solver_.set_preconditioner( parDGPreconditioner );
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
          solver_.unset_preconditioner();
        }
        else
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
      }

      unsigned int iterations () const
      {
        return solver_.number_of_iterations();
      }

    private:
      void setupSolver ( double redEps, double absLimit, unsigned int maxIterations, bool verbose )
      {
        PARDG::set_tolerance(solver_,redEps, absLimit, "fem.solver.errormeasure" );

        maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations );
        solver_.set_max_number_of_iterations( int( maxIterations ) );

        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose && Parameter :: verbose() )
        {
          solver_.IterativeSolver::set_output( std::cout );
          solver_.DynamicalObject::set_output( std::cout );
        }
      }

      mutable PARDG::BICGSTAB solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifdef USE_PARDG_ODE_SOLVER

#endif // #ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
