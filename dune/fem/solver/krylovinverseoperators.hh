#ifndef DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_INVERSEOPERATORS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/solver/cginverseoperator.hh>

#include <dune/fem/solver/linear/gmres.hh>
#include <dune/fem/solver/linear/bicgstab.hh>
#include <dune/fem/solver/linear/cg.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // KrylovInverseOperator
    // ---------------------

    template< class DiscreteFunction >
    class KrylovInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      typedef typename MPIManager::CollectiveCommunication  CollectiveCommunicationType;
    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op,
                                         double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op, double redEps, double absLimit,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, nullptr, redEps, absLimit,
                                          std::numeric_limits< unsigned int >::max(), readVerbose( parameter ), parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op, double redEps, double absLimit,
                                         unsigned int maxIterations,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, nullptr, redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                          readVerbose( parameter ), parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                         double redEps, double absLimit, unsigned int maxIterations,
                                         const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}

      KrylovInverseOperator ( double redEps, double absLimit,
                              unsigned int maxIterations, bool verbose,
                              const ParameterReader &parameter = Parameter::container() )
      : precondObj_(),
        tolerance_( absLimit ),
        errorType_( getErrorMeasure( parameter, "fem.solver.errormeasure" ) ),
        maxIterations_( std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations ) ),
        numOfIterations_( 0 ),
        verbose_( readVerbose( parameter, "fem.solver.verbose" ) ),
        method_( getMethod( parameter, "fem.solver.krylovmethod" ) ),
        restart_( method_ == gmres ? parameter.getValue< int >( "fem.solver.gmres.restart", 20 ) : 0 )
      {}

      KrylovInverseOperator ( double redEps, double absLimit,
                              const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( redEps, absLimit,
                               std::numeric_limits< unsigned int >::max(), readVerbose( parameter ), parameter ) {}

      KrylovInverseOperator ( double redEps, double absLimit,
                              unsigned int maxIterations,
                              const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( redEps, absLimit, maxIterations, readVerbose( parameter ), parameter ) {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        std::ostream* os = nullptr;
        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose_ && Parameter :: verbose() )
        {
          os = &std::cout;
        }


        int numIter = 0;

        if( method_ == gmres )
        {
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
          numIter = LinearSolver::gmres( *operator_, preconditioner_,
                                         v_, w, u, restart_,
                                         tolerance_, maxIterations_,
                                         errorType_, os );
        }
        else if( method_ == bicgstab )
        {
          if( v_.empty() )
          {
            v_.emplace_back( DomainFunctionType( "BiCGStab::r",   u.space() ) );
            v_.emplace_back( DomainFunctionType( "BiCGStab::r*",  u.space() ) );
            v_.emplace_back( DomainFunctionType( "BiCGStab::p",   u.space() ) );
            v_.emplace_back( DomainFunctionType( "BiCGStab::s",   u.space() ) );
            v_.emplace_back( DomainFunctionType( "BiCGStab::tmp", u.space() ) );
            if( preconditioner_ )
              v_.emplace_back( DomainFunctionType( "BiCGStab::z", u.space() ) );
          }

          // if solver convergence failed numIter will be negative
          numIter = LinearSolver::bicgstab( *operator_, preconditioner_,
                                            v_, w, u,
                                            tolerance_, maxIterations_,
                                            errorType_, os );
        }
        else if( method_ == cg )
        {
          if( v_.empty() )
          {
            v_.emplace_back( DomainFunctionType( "CG::h",   u.space() ) );
            v_.emplace_back( DomainFunctionType( "CG::r",  u.space() ) );
            v_.emplace_back( DomainFunctionType( "CG::p",   u.space() ) );

            if( preconditioner_ )
            {
              v_.emplace_back( DomainFunctionType( "CG::s",   u.space() ) );
              v_.emplace_back( DomainFunctionType( "CG::q", u.space() ) );
            }
          }

          // if solver convergence failed numIter will be negative
          numIter = LinearSolver::cg( *operator_, preconditioner_,
                                      v_, w, u,
                                      tolerance_, maxIterations_,
                                      errorType_, os );
        }

        // only add number of iterations when solver converged
        if( numIter > 0 )
          numOfIterations_ += numIter;
      }

      unsigned int iterations () const
      {
        return numOfIterations_;
      }

      void bind ( const OperatorType &op )
      {
        unbind();

        operator_ = &op;
        // if internal preconditioner is active set to preconditioner
        if( precondObj_ )
        {
          preconditioner_ = precondObj_.operator->();
        }
      }

      void bind ( const OperatorType &op, const PreconditionerType& preconditioner )
      {
        operator_ = &op;
        preconditioner_ = &preconditioner;
      }

      void unbind () { operator_ = nullptr; preconditioner_ = nullptr;  numOfIterations_ = 0; }

      void setMaxIterations ( unsigned int maxIterations )
      {
        maxIterations_ = std::min( static_cast< unsigned int >( std::numeric_limits< int >::max() ), maxIterations );
      }

    private:
      int getErrorMeasure( const ParameterReader& parameter, const char* paramName ) const
      {
        const std::string errorTypeTable[] =
          { "absolute", "relative", "residualreduction" };
        const int errorType = parameter.getEnum( paramName, errorTypeTable, 0 );
        return errorType ;
      }

      bool readVerbose( const ParameterReader& parameter, const char* paramName = "fem.solver.verbose" ) const
      {
        return parameter.getValue< bool >( paramName, false );
      }

      static const int cg       = 0 ; // CG
      static const int bicgstab = 1 ; // BiCGStab
      static const int gmres    = 2 ; // GMRES

      int getMethod( const ParameterReader& parameter, const char* paramName ) const
      {
        const std::string krylovMethodTable[] =
          { "cg", "bicgstab", "gmres" };
        const int methodType = parameter.getEnum( paramName, krylovMethodTable, gmres );
        return methodType;
      }

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op,
                              const PreconditionerType *preconditioner,
                              double redEps, double absLimit,
                              unsigned int maxIterations, bool verbose,
                              const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperator( redEps, absLimit, maxIterations, verbose, parameter )
      {
        bind(op);
        if( ! preconditioner_ )
        {
          const bool preconditioning = parameter.template getValue< bool >( "fem.preconditioning", false );
          if( preconditioning && std::is_base_of< AssembledOperator< DomainFunctionType, DomainFunctionType >, LinearOperator > :: value )
          {
            // create diagonal preconditioner
            precondObj_.reset( new DiagonalPreconditioner< DomainFunctionType, LinearOperator >( op ) );
            preconditioner_ = precondObj_.operator->();
          }
        }
      }

      const OperatorType *operator_ = nullptr;
      std::unique_ptr< PreconditionerType > precondObj_;
      const PreconditionerType *preconditioner_;

      mutable std::vector< DomainFunctionType > v_;

      const double tolerance_;
      const int errorType_;

      unsigned int maxIterations_;
      mutable unsigned int numOfIterations_;

      const bool verbose_;

      const int method_;
      const int restart_;
    };

  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
