#ifndef DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_INVERSEOPERATORS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/solver/cginverseoperator.hh>

#include <dune/fem/solver/linear/gmres.hh>
#include <dune/fem/solver/linear/bicgstab.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

#if 0
    // Internal Forward Declarations
    // -----------------------------

    namespace LinearSolver
    {
      template< class DomainFunction, class RangeFunction = DomainFunction >
      class OperatorAdapter;

      // OperatorAdapter for AdaptiveDiscreteFunction
      // ------------------------------------------

      template< class DomainFunctionSpace, class RangeFunctionSpace >
      class OperatorAdapter< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > >
      : public FunctionIF< typename RangeFunctionSpace :: RangeFieldType >
      {
        typedef OperatorAdapter< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > > ThisType;

        typedef AdaptiveDiscreteFunction< DomainFunctionSpace > DomainFunctionType;
        typedef AdaptiveDiscreteFunction< RangeFunctionSpace > RangeFunctionType;

        typedef typename RangeFunctionSpace :: RangeFieldType  RangeFieldType;

      public:
        typedef Operator< DomainFunctionType, RangeFunctionType > OperatorType;

        typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

        OperatorAdapter ( const OperatorType &op, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
        : operator_( op ),
          domainSpace_( domainSpace ),
          rangeSpace_( rangeSpace )
        {}

        void operator() ( const RangeFieldType *u, RangeFieldType *w )
        {
          DomainFunctionType uFunction( "OperatorAdapter u", domainSpace_, (const RangeFieldType *) u );
          RangeFunctionType  wFunction( "OperatorAdapter w", rangeSpace_ , (RangeFieldType *) w );
          operator_( uFunction, wFunction );
        }

        int size() const
        {
          return domainSpace_.size();
        }

      private:
        const OperatorType &operator_;
        const DomainFunctionSpaceType &domainSpace_;
        const RangeFunctionSpaceType &rangeSpace_;
      };
    }  // end namespace LinearSolver

#endif

    // GeneralizedMinResInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    class GeneralizedMinResInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      //typedef LinearSolver::OperatorAdapter< DiscreteFunction, DiscreteFunction > OperatorAdapterType;

      typedef typename MPIManager::CollectiveCommunication  CollectiveCommunicationType;
    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      GeneralizedMinResInverseOperator ( const OperatorType &op,
                                              double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                              const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      GeneralizedMinResInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                         const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit,
                                          std::numeric_limits< unsigned int >::max(), readVerbose( parameter ), parameter )
      {}

      GeneralizedMinResInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                         unsigned int maxIterations,
                                         const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}


      GeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                         const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      GeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit,
                                         const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                          readVerbose( parameter ), parameter )
      {}

      GeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit, unsigned int maxIterations,
                                         const ParameterReader &parameter = Parameter::container() )
      : GeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}


      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        std::ostream* os = nullptr;
        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose_ && Parameter :: verbose() )
        {
          os = &std::cout;
        }

        // TODO: Allow to store memory in class, this either needs pointer or
        // operator to export the space
        if( v_.empty() )
        {
          v_.reserve( restart_+1 );
          for( int i=0; i<=restart_; ++i )
          {
            v_.emplace_back( DiscreteFunction( "GMRes::v", u.space() ) );
          }
          if( preconditioner_ )
            v_.emplace_back( DiscreteFunction( "GMRes::z", u.space() ) );
        }

        // if solver convergence failed numIter will be negative
        int numIter = LinearSolver::gmres( operator_, preconditioner_,
                                           v_, w, u, restart_,
                                           tolerance_, maxIterations_,
                                           errorType_, os );

        // only add number of iterations when solver converged
        if( numIter > 0 )
          numOfIterations_ += numIter;
      }

      unsigned int iterations () const
      {
        return numOfIterations_;
      }

    private:
      bool readVerbose( const ParameterReader& parameter ) const
      {
        return parameter.getValue< bool >("fem.solver.verbose", false );
      }

      GeneralizedMinResInverseOperator ( const OperatorType &op,
                                         const PreconditionerType *preconditioner,
                                         double redEps, double absLimit,
                                         unsigned int maxIterations, bool verbose,
                                         const ParameterReader &parameter = Parameter::container() )
      : //  solver_( MPIManager::comm(), parameter.getValue< int >( "fem.solver.gmres.restart", 20 ) ),
        operator_( op ),
        preconditioner_( preconditioner ),
        tolerance_( absLimit ),
        errorType_( LinearSolver::getErrorMeasure( parameter, "fem.solver.errormeasure" ) ),
        maxIterations_( std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations ) ),
        numOfIterations_( 0 ),
        verbose_( parameter.getValue< bool >("fem.solver.verbose", false ) ),
        restart_( parameter.getValue< int >( "fem.solver.gmres.restart", 20 ) )
      {
      }

      // mutable LinearSolver::GMRES<double, CollectiveCommunicationType > solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;

      mutable std::vector< DomainFunctionType > v_;

      const double tolerance_;
      const int errorType_;

      const unsigned int maxIterations_;
      mutable unsigned int numOfIterations_;

      const bool verbose_;

      const int restart_;
    };

    // BiCGStabInverseOperator
    // ----------------------------

    template< class DiscreteFunction >
    class BiCGStabInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      // typedef LinearSolver::OperatorAdapter< DiscreteFunction, DiscreteFunction > OperatorAdapterType;

    public:
      typedef typename MPIManager::CollectiveCommunication  CollectiveCommunicationType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      BiCGStabInverseOperator ( const OperatorType &op,
                                double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      BiCGStabInverseOperator ( const OperatorType &op,
                                double redEps, double absLimit,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, nullptr, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                 readVerbose( parameter ), parameter ) {}

      BiCGStabInverseOperator ( const OperatorType &op,
                                double redEps, double absLimit, unsigned int maxIterations,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, nullptr, redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}


      BiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      BiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                double redEps, double absLimit,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                 readVerbose( parameter ), parameter ) {}

      BiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                double redEps, double absLimit, unsigned int maxIterations,
                                const ParameterReader &parameter = Parameter::container() )
      : BiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, false, parameter ) {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        std::ostream* os = nullptr;
        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose_ && Parameter :: verbose() )
        {
          os = &std::cout;
        }

        // TODO: Allow to store memory in class, this either needs pointer or
        // operator to export the space

        if( tempMemory_.empty() )
        {
          tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::r",   u.space() ) );
          tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::r*",  u.space() ) );
          tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::p",   u.space() ) );
          tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::s",   u.space() ) );
          tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::tmp", u.space() ) );
          if( preconditioner_ )
            tempMemory_.emplace_back( DomainFunctionType( "BiCGStab::z", u.space() ) );
        }

        // if solver convergence failed numIter will be negative
        int numIter = LinearSolver::bicgstab( operator_, preconditioner_,
                                              tempMemory_, w, u,
                                              tolerance_, maxIterations_,
                                              errorType_, os );

        // only add number of iterations when solver converged
        if( numIter > 0 )
          numOfIterations_ += numIter;
      }

      unsigned int iterations () const
      {
        return numOfIterations_;
      }

    private:
      bool readVerbose( const ParameterReader &parameter ) const
      {
        return parameter.getValue< bool >("fem.solver.verbose", false );
      }

      BiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType *preconditioner,
                                double redEps, double absLimit, unsigned int maxIterations, bool verb,
                                const ParameterReader &parameter )
      : //solver_( MPIManager::comm() ),
        operator_( op ),
        preconditioner_( preconditioner ),
        tolerance_( absLimit ),
        errorType_( LinearSolver::getErrorMeasure( parameter, "fem.solver.errormeasure" ) ),
        maxIterations_( std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations ) ),
        numOfIterations_( 0 ),
        verbose_( parameter.getValue< bool >("fem.solver.verbose", false ) )
      {
      }

      // mutable LinearSolver::BICGSTAB<double, CollectiveCommunicationType > solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;

      mutable std::vector< DomainFunctionType > tempMemory_;

      const double tolerance_;
      const int errorType_;

      const unsigned int maxIterations_;
      mutable unsigned int numOfIterations_;

      const bool verbose_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
