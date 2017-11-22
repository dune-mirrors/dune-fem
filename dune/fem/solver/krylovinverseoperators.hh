#ifndef DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_INVERSEOPERATORS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/parameter.hh>

#include <dune/fem/solver/cginverseoperator.hh>

#include <dune/fem/solver/linear/gmres.hh>
#include <dune/fem/solver/linear/bicgstab.hh>
#include <dune/fem/solver/linear/cg.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // KrylovInverseOperatorBase
    // ---------------------

    template< class DiscreteFunction, int method = -1 >
    class KrylovInverseOperatorBase
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
      KrylovInverseOperatorBase ( const LinearOperator &op,
                                  double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op, double redEps, double absLimit,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, nullptr, redEps, absLimit,
                                   parameter.maxLinearIterationsParameter(), parameter.verbose() ) {}

      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op, double redEps, double absLimit,
                                  unsigned int maxIterations,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, nullptr, redEps, absLimit, maxIterations, parameter.verbose() ) {}

      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                  double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                  double redEps, double absLimit,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, &preconditioner, redEps, absLimit, parameter.maxLinearIterationsParameter(),
                                   parameter.verbose() ) {}

      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op, const PreconditionerType &preconditioner,
                                  double redEps, double absLimit, unsigned int maxIterations,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( op, &preconditioner, redEps, absLimit, maxIterations, parameter.verbose() ) {}

      KrylovInverseOperatorBase ( double redEps, double absLimit,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( redEps, absLimit, parameter.maxLinearIterationsParameter(), parameter.verbose() ) {}

      KrylovInverseOperatorBase ( double redEps, double absLimit,
                                  unsigned int maxIterations,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperatorBase( redEps, absLimit, maxIterations,  parameter.verbose() ) {}

      KrylovInverseOperatorBase ( double redEps, double absLimit,
                                  unsigned int maxIterations, bool verbose,
                                  const ParameterReader& parameter )
        : KrylovInverseOperatorBase( redEps, absLimit, maxIterations, verbose,
            SolverParameter( parameter ) ) {}

      //! main constructor
      KrylovInverseOperatorBase ( double redEps, double absLimit,
                                  unsigned int maxIterations, bool verbose,
                                  const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : precondObj_(),
        tolerance_( absLimit ),
        errorType_( parameter.errorMeasure() ),
        maxIterations_( std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations ) ),
        numOfIterations_( 0 ),
        verbose_( verbose ? true : parameter.verbose() ), // verbose overrules parameter.verbose()
        method_( method < 0 ? parameter.krylovMethod() : method ),
        restart_( method_ == SolverParameter::gmres ? parameter.gmresRestart() : 0 )
      {}

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

        if( method_ == SolverParameter::gmres )
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
        else if( method_ == SolverParameter::bicgstab )
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
        else if( method_ == SolverParameter::cg )
        {
          if( v_.empty() )
          {
            v_.emplace_back( DomainFunctionType( "CG::h",  u.space() ) );
            v_.emplace_back( DomainFunctionType( "CG::r",  u.space() ) );
            v_.emplace_back( DomainFunctionType( "CG::p",  u.space() ) );

            if( preconditioner_ )
            {
              v_.emplace_back( DomainFunctionType( "CG::s", u.space() ) );
              v_.emplace_back( DomainFunctionType( "CG::q", u.space() ) );
            }
          }

          // if solver convergence failed numIter will be negative
          numIter = LinearSolver::cg( *operator_, preconditioner_,
                                      v_, w, u,
                                      tolerance_, maxIterations_,
                                      errorType_, os );
        }

        // only store number of iterations when solver converged, otherwise numIter < 0
        numOfIterations_ = ( numIter > 0 ) ? numIter : 0;
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
      template <class LinearOperator>
      KrylovInverseOperatorBase ( const LinearOperator &op,
                              const PreconditionerType *preconditioner,
                              double redEps, double absLimit,
                              unsigned int maxIterations, bool verbose,
                              const ParameterReader &parameter = Parameter::container() )
      : KrylovInverseOperatorBase( redEps, absLimit, maxIterations, verbose, parameter )
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


    // KrylovInverseOperator
    // ---------------------

    template< class DiscreteFunction >
    using KrylovInverseOperator = KrylovInverseOperatorBase< DiscreteFunction, -1 >;


    // CgInverseOperator
    // -----------------

    template< class DiscreteFunction >
    using CgInverseOperator = KrylovInverseOperatorBase< DiscreteFunction, SolverParameter :: cg >;


    // BicgstabInverseOperator
    // -----------------------

    template< class DiscreteFunction >
    using BicgstabInverseOperator = KrylovInverseOperatorBase< DiscreteFunction, SolverParameter :: bicgstab >;


    // GmresInverseOperator
    // --------------------

    template< class DiscreteFunction >
    using GmresInverseOperator = KrylovInverseOperatorBase< DiscreteFunction, SolverParameter :: gmres >;


    // ParDGGeneralizedMinResInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    using ParDGGeneralizedMinResInverseOperator = GmresInverseOperator< DiscreteFunction >;

    // ParDGBiCGStabInverseOperator
    // ----------------------------

    template< class DiscreteFunction >
    using ParDGBiCGStabInverseOperator = BicgstabInverseOperator< DiscreteFunction >;


  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
