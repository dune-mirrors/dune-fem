#ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/solver/parameter.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>

#if HAVE_DUNE_ISTL
#include <dune/common/version.hh>

#include <dune/fem/operator/linear/istladapter.hh>
#include <dune/fem/operator/linear/istloperator.hh>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioner.hh>

namespace Dune
{

  namespace Fem
  {

    // wrapper for Fem Preconditioners (Operators acting as preconditioners) into ISTL preconditioners
    template< class Preconditioner >
    class ISTLPreconditionAdapter
     : public Dune::Preconditioner< typename Preconditioner::RangeFunctionType::DofStorageType, typename Preconditioner::DomainFunctionType::DofStorageType >
    {
      typedef ISTLPreconditionAdapter< Preconditioner > ThisType;
      typedef Dune::Preconditioner< typename Preconditioner::RangeFunctionType::DofStorageType, typename Preconditioner::DomainFunctionType::DofStorageType > BaseType;

      typedef typename Preconditioner::DomainFunctionType DomainFunctionType;
      typedef typename Preconditioner::RangeFunctionType  RangeFunctionType;

    public:
      typedef typename BaseType::domain_type domain_type;
      typedef typename BaseType::range_type  range_type;
      typedef typename BaseType::field_type  field_type;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType  RangeFunctionSpaceType;

      ISTLPreconditionAdapter ( const Preconditioner *precon, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
      : precon_( precon ),
        domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace )
      {}

      // pre and post do nothing here
      virtual void pre ( domain_type &x, range_type &y ) override {}
      virtual void post ( domain_type &x ) override {}

      virtual void apply ( domain_type &x, const range_type &y ) override
      {
        // no precon
        if( !precon_ )
        {
          x = y;
        }
        else
        {
          // note: ISTL switches the arguments !!!
          // it is assumed that we have a left preconditioner
          RangeFunctionType px( "ISTLPreconditionAdapter::apply::x", rangeSpace_, x );
          DomainFunctionType py( "ISTLPreconditionAdapter::apply::y", domainSpace_, y );

          (*precon_)( px, py );
        }
      }

      SolverCategory::Category category () const { return SolverCategory::sequential; }

    protected:
      const Preconditioner *precon_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };


    template< class BlockVector >
    struct ISTLSolverReduction
    {
      ISTLSolverReduction ( std::shared_ptr< ISTLSolverParameter > parameter )
        : parameter_(  parameter )
      {}

      double operator() ( const Dune::LinearOperator< BlockVector, BlockVector > &op,
                          Dune::ScalarProduct< BlockVector > &scp,
                          const BlockVector &rhs, const BlockVector &x ) const
      {

        if( parameter_->errorMeasure() == 0)
        {
          BlockVector residuum( rhs );
          op.applyscaleadd( -1., x, residuum );
          const double res = scp.norm( residuum );
          return (res > 0 ? parameter_->tolerance() / res : 1e-3);
        }
        else
          return parameter_->tolerance();
      }

    private:
      std::shared_ptr<ISTLSolverParameter> parameter_;
    };


    struct ISTLInverseOperatorMethods
    {
      static std::vector< int > supportedSolverMethods() {
        return std::vector< int > ({
                                     SolverParameter::gmres, // default solver
                                     SolverParameter::cg,
                                     SolverParameter::bicgstab,
                                     SolverParameter::minres,
                                     SolverParameter::gradient,
                                     SolverParameter::loop,
                                     SolverParameter::superlu
                                   });
      }
    };

    template< int method,
              class X,
              class Reduction = ISTLSolverReduction< X > >
    struct ISTLSolverAdapter
    {
      typedef Reduction ReductionType;

      typedef X domain_type;
      typedef X range_type;

      ISTLSolverAdapter ( const ReductionType &reduction, std::shared_ptr<ISTLSolverParameter> parameter )
        : reduction_( reduction ),
          method_( method < 0 ? parameter->solverMethod( ISTLInverseOperatorMethods::supportedSolverMethods() ) : method ),
          parameter_( parameter )
      {}

      template<class Op, class ScP, class PC >
      void operator () ( Op& op, ScP &scp, PC &pc,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
        int verbosity = (Parameter::verbose( Parameter::solverStatistics ) && parameter_->verbose()) ? 2 : 0;
        if ( verbosity && Parameter::verbose( Parameter::extendedStatistics ) )
          verbosity += 2;

        int maxIterations = std::min( std::numeric_limits< int >::max(), parameter_->maxIterations() );
        if( method_ == SolverParameter::cg )
        {
          typedef Dune::CGSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::bicgstab )
        {
          typedef Dune::BiCGSTABSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::gmres )
        {
          typedef Dune::RestartedGMResSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), parameter_->gmresRestart(), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::minres )
        {
          typedef Dune::MINRESSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
        }
        else if( method_ == SolverParameter::gradient )
        {
          typedef Dune::GradientSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::loop )
        {
          typedef Dune::LoopSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::superlu )
        {
          callSuperLU( op, rhs, x, result );
        }
        else
        {
          DUNE_THROW(NotImplemented,"ISTLSolverAdapter::operator(): wrong method solver identifier" << method_ );
        }
      }

      void setMaxLinearIterations( unsigned int maxIterations ) { parameter_->setMaxIterations(maxIterations); }
      void setMaxIterations( unsigned int maxIterations ) { parameter_->setMaxIterations(maxIterations); }

      std::shared_ptr<ISTLSolverParameter> parameter () const { return parameter_; }

    protected:
      template< class ImprovedMatrix >
      void callSuperLU ( ISTLParallelMatrixAdapterInterface< ImprovedMatrix >& op,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
#if HAVE_SUPERLU
        const int verbosity = (Dune::Fem::Parameter::verbose( Parameter::solverStatistics ) && parameter_->verbose()) ? 2 : 0;
        typedef typename ImprovedMatrix :: BaseType Matrix;
        const ImprovedMatrix& matrix = op.getmat();
        SuperLU< Matrix > solver( matrix, verbosity );
        solver.apply( x, rhs, result );
#else
        DUNE_THROW(NotImplemented,"ISTLSolverAdapter::callSuperLU: SuperLU solver selected but SuperLU not available!");
#endif
      }

      template< class Op >
      void callSuperLU ( ISTLLinearOperatorAdapter< Op >& op, range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
        DUNE_THROW(NotImplemented,"ISTLSolverAdapter::callSuperLU: SuperLU only works for AssembledLinearOperators!");
      }

      ReductionType reduction_;
      const int method_;

      std::shared_ptr<ISTLSolverParameter> parameter_;
    };



    // ISTLInverseOperator
    // -------------------

    template< class DiscreteFunction, int method = -1,
              class Preconditioner = Fem::Operator< DiscreteFunction, DiscreteFunction > >
    class ISTLInverseOperator;

    template< class DiscreteFunction, int method, class Preconditioner >
    struct ISTLInverseOperatorTraits
    {
      typedef DiscreteFunction  DiscreteFunctionType;
      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Preconditioner PreconditionerType;

      typedef ISTLBlockVectorDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType ;
      typedef Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;

      typedef ISTLInverseOperator< DiscreteFunction, method, Preconditioner >  InverseOperatorType;
      typedef ISTLSolverParameter SolverParameterType;
    };

    template< class DiscreteFunction, int method, class Preconditioner >
    class ISTLInverseOperator : public InverseOperatorInterface<
                                ISTLInverseOperatorTraits< DiscreteFunction,
                                method, Preconditioner > >
    {
      typedef ISTLInverseOperatorTraits< DiscreteFunction, method, Preconditioner > Traits;
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
      using BaseType :: setMaxIterations;
      using BaseType :: setMaxLinearIterations;

    protected:
      typedef typename Traits::SolverDiscreteFunctionType       SolverDiscreteFunctionType;

      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;


      typedef typename SolverDiscreteFunctionType::ScalarProductType   ParallelScalarProductType;
      typedef typename SolverDiscreteFunctionType::DofStorageType      BlockVectorType;

      typedef ISTLSolverAdapter< method, BlockVectorType > SolverAdapterType;
      typedef typename SolverAdapterType::ReductionType ReductionType;
    public:

      //non-deprecated constructors
      ISTLInverseOperator ( const ISTLSolverParameter & parameter = ISTLSolverParameter() )
        : BaseType( parameter ), solverAdapter_( ReductionType( parameter_ ), parameter_ )
      {}

      ISTLInverseOperator ( const OperatorType &op,
                            const ISTLSolverParameter & parameter = ISTLSolverParameter() )
        : ISTLInverseOperator ( parameter )
      {
        bind( op );
      }

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            const ISTLSolverParameter & parameter = ISTLSolverParameter() )
        : ISTLInverseOperator( parameter )
      {
        bind( op, preconditioner );
      }

    protected:
      // apply for arbitrary domain function type and matching range function type
      template <class DomainFunction>
      int apply( const DomainFunction& u, SolverDiscreteFunctionType& w ) const
      {
        auto& scp = w.scalarProduct();
        // u may not be a discrete function, therefore use w.space()
        const DiscreteFunctionSpaceType& space = w.space();

        typedef Dune::Preconditioner< BlockVectorType, BlockVectorType > ISTLPreconditionerType;
        std::shared_ptr< ISTLPreconditionerType > istlPre;
        if constexpr (std::is_same< SolverDiscreteFunctionType, RangeFunctionType >::value )
        {
          typedef ISTLPreconditionAdapter< PreconditionerType > ISTLPreconditionerAdapterType;
          istlPre.reset( new ISTLPreconditionerAdapterType( preconditioner_, space, space ) );
        }

        if( assembledOperator_ )
        {
          auto& matrix = assembledOperator_->matrixAdapter( *(solverAdapter_.parameter()) );
          // if preconditioner_ was set use that one, otherwise the one from the matrix object
          ISTLPreconditionerType& matrixPre = matrix.preconditionAdapter();
          ISTLPreconditionerType& precon    = ( preconditioner_ ) ? (*istlPre) : matrixPre;
          return solve( matrix, scp, precon, u, w );
        }

        if constexpr (std::is_same< SolverDiscreteFunctionType, RangeFunctionType >::value )
        {
          assert( operator_ );
          assert( istlPre );
          typedef ISTLLinearOperatorAdapter< OperatorType >  ISTLOperatorType;
          ISTLOperatorType istlOperator( *operator_, space, space );
          return solve( istlOperator, scp, *istlPre, u, w );
        }

        DUNE_THROW(InvalidStateException,"ISTLInverseOperator::apply: No valid operator found!");
        return -1;
      }

      //! final solve execution only copying the right hand side
      template< class OperatorAdapter, class ISTLPreconditioner, class DomainFunction >
      int solve ( OperatorAdapter &istlOperator, ParallelScalarProductType &scp,
                  ISTLPreconditioner &preconditioner,
                  const DomainFunction& u,
                  SolverDiscreteFunctionType& w ) const
      {
        if( ! rhs_ )
        {
          // u may not be a discrete function, therefore use w.space()
          rhs_.reset( new SolverDiscreteFunctionType( "ISTLInvOp::rhs", w.space() ) );
          rightHandSideCopied_ = false;
        }

        if( ! rightHandSideCopied_ )
        {
          // copy right hand side since ISTL solvers seem to modify it
          rhs_->assign( u );
          rightHandSideCopied_ = true;
        }

        Dune::InverseOperatorResult result;
        solverAdapter_.setMaxIterations( parameter_->maxIterations() );
        solverAdapter_( istlOperator, scp, preconditioner, rhs_->blockVector(), w.blockVector(), result );
        return (result.converged) ? result.iterations : -(result.iterations);
      }

      using BaseType :: operator_;
      using BaseType :: assembledOperator_;
      using BaseType :: preconditioner_;

      using BaseType :: rhs_;
      using BaseType :: x_;

      using BaseType :: rightHandSideCopied_;
      using BaseType :: parameter_;

      mutable SolverAdapterType solverAdapter_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
