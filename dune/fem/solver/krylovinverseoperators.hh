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
#include <dune/fem/solver/inverseoperatorinterface.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // KrylovInverseOperator
    // ---------------------

    template< class DiscreteFunction, int method = -1 >
    class KrylovInverseOperator;

    template< class DiscreteFunction, int method >
    struct KrylovInverseOperatorTraits
    {
      typedef DiscreteFunction    DiscreteFunctionType;
      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef OperatorType  PreconditionerType;

      typedef OperatorType AssembledOperatorType;
      typedef DiscreteFunction SolverDiscreteFunctionType ;

      typedef KrylovInverseOperator< DiscreteFunction, method >  InverseOperatorType;
    };


    template< class DiscreteFunction, int method >
    class KrylovInverseOperator
    : public InverseOperatorInterface< KrylovInverseOperatorTraits< DiscreteFunction, method > >
    {
      typedef KrylovInverseOperatorTraits< DiscreteFunction, method > Traits;
      typedef InverseOperatorInterface< Traits > BaseType;

      friend class InverseOperatorInterface< Traits >;
    public:
      typedef typename BaseType::DomainFunctionType     DomainFunctionType;
      typedef typename BaseType::RangeFunctionType      RangeFunctionType;
      typedef typename BaseType::OperatorType           OperatorType;
      typedef typename BaseType::PreconditionerType     PreconditionerType;
      typedef typename BaseType::AssembledOperatorType  AssembledOperatorType;

      using BaseType :: bind;
      using BaseType :: unbind;
      using BaseType :: setMaxLinearIterations;
      using BaseType :: setMaxIterations;

    public:
      //deprecated constructors
      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op, double redEps, double absLimit,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( op, redEps, absLimit,
                                   parameter.maxLinearIterationsParameter(), parameter.verbose(), parameter ) {}

      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op, double redEps, double absLimit,
                              unsigned int maxIterations,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( op, redEps, absLimit, maxIterations, parameter.verbose(), parameter ) {}

      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                              double redEps, double absLimit,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( op, &preconditioner, redEps, absLimit, parameter.maxLinearIterationsParameter(), parameter.verbose(), parameter) {}

      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                              double redEps, double absLimit, unsigned int maxIterations,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, parameter.verbose(), parameter ) {}

      [[deprecated]]
      KrylovInverseOperator ( double redEps, double absLimit,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( redEps, absLimit, parameter.maxLinearIterationsParameter(), parameter.verbose(), parameter ) {}


      [[deprecated]]
      KrylovInverseOperator ( double redEps, double absLimit,
                              unsigned int maxIterations,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( redEps, absLimit, maxIterations,  parameter.verbose(), parameter ) {}

      [[deprecated]]
      KrylovInverseOperator (  double redEps, double absLimit,
                               unsigned int maxIterations, bool verbose,
                               const SolverParameter &parameter = SolverParameter(Parameter::container()) )
        : KrylovInverseOperator( parameter )
      {
        parameter_->setReductionTol( redEps );
        parameter_->setAbsoluteTol( absLimit );
        parameter_->setMaxIterations( maxIterations );
        parameter_->setVerbose( verbose );
      }

      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op,
                              double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( parameter )
      {
        bind( op );
        parameter_->setReductionTol( redEps );
        parameter_->setAbsoluteTol( absLimit );
        parameter_->setMaxIterations( maxIterations );
        parameter_->setVerbose( verbose );
     }

      template <class LinearOperator>
      [[deprecated]]
      KrylovInverseOperator ( const LinearOperator &op, const PreconditionerType &preconditioner,
                              double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                              const SolverParameter &parameter = SolverParameter() )
      : KrylovInverseOperator( op, &preconditioner, parameter )
      {
        parameter_->setReductionTol( redEps );
        parameter_->setAbsoluteTol( absLimit );
        parameter_->setMaxIterations( maxIterations );
        parameter_->setVerbose( verbose );
      }

      [[deprecated]]
      KrylovInverseOperator ( double redEps, double absLimit,
                              unsigned int maxIterations, bool verbose,
                              const ParameterReader& parameter )
        : KrylovInverseOperator( redEps, absLimit, maxIterations, verbose,
            SolverParameter( parameter ) ) {}

      //non-deprecated constructors
      KrylovInverseOperator ( const OperatorType &op, const SolverParameter &parameter = SolverParameter(Parameter::container()) )
        : KrylovInverseOperator( parameter )
      {
        bind( op );
      }

      KrylovInverseOperator ( const OperatorType &op, const PreconditionerType& preconditioner,
                              const SolverParameter &parameter = SolverParameter(Parameter::container()) )
        : KrylovInverseOperator( op, &preconditioner, parameter )
      {}

      //! main constructor
      KrylovInverseOperator ( const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter ),
        precondObj_(),
        verbose_( parameter.verbose() ),
        method_( method < 0 ? parameter.krylovMethod() : method )
      {}

      void bind ( const OperatorType &op )
      {
        if( precondObj_ )
          BaseType::bind( op, *precondObj_ );
        else
          BaseType::bind( op );
      }

    protected:
      int apply( const DomainFunctionType &u, RangeFunctionType &w ) const
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
            v_.reserve( parameter_->gmresRestart()+1 );
            for( int i=0; i<=parameter_->gmresRestart(); ++i )
            {
              v_.emplace_back( DiscreteFunction( "GMRes::v", u.space() ) );
            }
            if( preconditioner_ )
              v_.emplace_back( DiscreteFunction( "GMRes::z", u.space() ) );
          }

          // if solver convergence failed numIter will be negative
          numIter = LinearSolver::gmres( *operator_, preconditioner_,
                                         v_, w, u, parameter_->gmresRestart(),
                                         parameter_->absoluteTol(), parameter_->maxIterations(),
                                         parameter_->errorMeasure(), os );
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
                                            parameter_->absoluteTol(), parameter_->maxIterations(),
                                            parameter_->errorMeasure(), os );
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
                                      parameter_->absoluteTol(), parameter_->maxIterations(),
                                      parameter_->errorMeasure(), os );
        }
        else
        {
          DUNE_THROW(InvalidStateException,"KrylovInverseOperator: invalid method " << method_ );
        }

        return numIter;
      }

      template <class LinearOperator>
      KrylovInverseOperator ( const LinearOperator &op,
                              const PreconditionerType *preconditioner,
                              const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : KrylovInverseOperator( parameter )
      {
        bind(op);
        if( ! preconditioner_ )
        {
          const bool preconditioning = parameter.parameter().template getValue< bool >( "fem.preconditioning", false );
          if( preconditioning && std::is_base_of< AssembledOperator< DomainFunctionType, DomainFunctionType >, LinearOperator > :: value )
          {
            // create diagonal preconditioner
            precondObj_.reset( new DiagonalPreconditioner< DomainFunctionType, LinearOperator >( op ) );
            preconditioner_ = precondObj_.operator->();
          }
        }
      }

    protected:
      using BaseType :: operator_;
      using BaseType :: preconditioner_;

      using BaseType :: parameter_;
      using BaseType :: iterations_;

      std::unique_ptr< PreconditionerType > precondObj_;

      mutable std::vector< DomainFunctionType > v_;

      const bool verbose_;

      const int method_;
    };


    // CgInverseOperator
    // -----------------

    template< class DiscreteFunction >
    using CgInverseOperator = KrylovInverseOperator< DiscreteFunction, SolverParameter :: cg >;


    // BicgstabInverseOperator
    // -----------------------

    template< class DiscreteFunction >
    using BicgstabInverseOperator = KrylovInverseOperator< DiscreteFunction, SolverParameter :: bicgstab >;


    // GmresInverseOperator
    // --------------------

    template< class DiscreteFunction >
    using GmresInverseOperator = KrylovInverseOperator< DiscreteFunction, SolverParameter :: gmres >;


    // ParDGGeneralizedMinResInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    using ParDGGeneralizedMinResInverseOperator = GmresInverseOperator< DiscreteFunction >;

    // ParDGBiCGStabInverseOperator
    // ----------------------------

    template< class DiscreteFunction >
    using ParDGBiCGStabInverseOperator = BicgstabInverseOperator< DiscreteFunction >;

    template <class DiscreteFunctionType, class OpType >
    using OEMCGOp = CgInverseOperator< DiscreteFunctionType >;

    template <class DiscreteFunctionType, class OpType >
    using OEMBICGSTABOp = BicgstabInverseOperator< DiscreteFunctionType >;

    template <class DiscreteFunctionType, class OpType >
    using OEMBICGSQOp = BicgstabInverseOperator< DiscreteFunctionType >;

    template <class DiscreteFunctionType, class OpType >
    using OEMGMRESOp = GmresInverseOperator< DiscreteFunctionType >;

    template <class DiscreteFunctionType, class OpType >
    using GMRESOp = GmresInverseOperator< DiscreteFunctionType >;

  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_SOLVER_INVERSEOPERATORS_HH
