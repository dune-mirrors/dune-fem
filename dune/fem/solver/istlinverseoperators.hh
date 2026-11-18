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

#if HAVE_DUNE_ISTL_GPU
#include <dune/istl/gpu/vector.hh>
#include <dune/istl/gpu/sparsematrix.hh>
#include <dune/istl/gpu/hyprepreconditioner.hh>
#endif

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

      SolverCategory::Category category () const override { return SolverCategory::sequential; }

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
                                     SolverParameter::superlu,
                                     SolverParameter::fgmres
                                   });
      }
    };

    template< int method,
              class X,
              class Xgpu,
              class Reduction = ISTLSolverReduction< Xgpu > >
    struct ISTLSolverAdapter
    {
      typedef Reduction ReductionType;

      typedef Xgpu domain_type;
      typedef Xgpu range_type;

      // if X and Xgpu differ a copy is needed
      static const bool copyVector = ! std::is_same_v< X, Xgpu >;

      ISTLSolverAdapter ( const ReductionType &reduction, std::shared_ptr<ISTLSolverParameter> parameter )
        : reduction_( reduction ),
          method_( method < 0 ? parameter->solverMethod( ISTLInverseOperatorMethods::supportedSolverMethods() ) : method ),
          parameter_( parameter )
      {}

      template<class Op, class ScP, class PC, class VectorType >
      void operator () ( Op& op, ScP &scp, PC &pc,
                         VectorType &rhs, VectorType &x,
                         Dune::InverseOperatorResult &result ) const
      {
        if constexpr ( std::is_same_v< VectorType, domain_type > )
        {
          applySolver( op, scp, pc, rhs, x, result );
        }
        else
        {
          std::cout << "Copying to device memory" << std::endl;
          // otherwise copy to device, solve, copy back
          const size_t size = rhs.size()*X::block_type::dimension ;
          range_type gpuRhs( size );
          gpuRhs.copyFromHost( rhs );
          domain_type gpuX( size );
          gpuX.copyFromHost( x );

          typedef Dune::GPU::GpuSparseMatrix< typename domain_type::field_type > GpuSparseMatrixType;
          GpuSparseMatrixType gpuMat = GpuSparseMatrixType::fromMatrix( op.getmat() /*, true */ );
          using GPUOperator = Dune::MatrixAdapter< GpuSparseMatrixType, domain_type, domain_type >;
          GPUOperator gpuOp( gpuMat );
          Dune::ScalarProduct< domain_type > gpuScp;
          Dune::Richardson< domain_type, domain_type > gpuPc;
          applySolver( gpuOp, gpuScp, gpuPc, gpuRhs, gpuX, result );
          gpuX.copyToHost( x );
        }
      }

      template<class Op, class ScP, class PC, class VectorType >
      void applySolver(Op& op, ScP &scp, PC &pc,
                       VectorType &rhs, VectorType &x,
                       Dune::InverseOperatorResult &result ) const
      {
        int verbosity = (Parameter::verbose( Parameter::solverStatistics ) && parameter_->verbose()) ? 2 : 0;
        if ( verbosity && Parameter::verbose( Parameter::extendedStatistics ) )
          verbosity += 2;

        int maxIterations = std::min( std::numeric_limits< int >::max(), parameter_->maxIterations() );
        if( method_ == SolverParameter::cg )
        {
          typedef Dune::CGSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::bicgstab )
        {
          typedef Dune::BiCGSTABSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::gmres )
        {
          typedef Dune::RestartedGMResSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), parameter_->gmresRestart(), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::fgmres )
        {
          typedef Dune::RestartedFlexibleGMResSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), parameter_->gmresRestart(), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::minres )
        {
          typedef Dune::MINRESSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
        }
        else if( method_ == SolverParameter::gradient )
        {
          typedef Dune::GradientSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::loop )
        {
          typedef Dune::LoopSolver< VectorType > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbosity );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::superlu )
        {
#ifndef USE_IBCRS_MATRIX
          callSuperLU( op, rhs, x, result );
#else
          DUNE_THROW(NotImplemented,"SuperLU not supported for Experimental::BCRSMatrix");
#endif
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
      void callSuperLU ( Op& op, range_type &rhs, domain_type &x,
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

      typedef ISTLSolverAdapter< method, BlockVectorType, BlockVectorType > SolverAdapterType;

#if HAVE_DUNE_ISTL_GPU
      typedef Dune::GPU::GpuVector< typename BlockVectorType::field_type > GPUVectorType;
      typedef ISTLSolverAdapter< method, BlockVectorType, GPUVectorType > GPUSolverAdapterType;
#else
      typedef SolverAdapterType GPUSolverAdapterType;
#endif

    public:

      //non-deprecated constructors
      ISTLInverseOperator ( const ISTLSolverParameter & parameter = ISTLSolverParameter() )
        : BaseType( parameter ),
          solverAdapter_(),
          gpuSolverAdapter_()
      {
#if HAVE_DUNE_ISTL_GPU
        if( parameter.usegpu() )
        {
          typedef typename GPUSolverAdapterType::ReductionType ReductionType;
          gpuSolverAdapter_.reset( new GPUSolverAdapterType( ReductionType( parameter_ ), parameter_  ));
        }
#endif

        {
          typedef typename SolverAdapterType::ReductionType ReductionType;
          solverAdapter_.reset( new SolverAdapterType( ReductionType( parameter_ ), parameter_ ) );
        }
      }

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
#if HAVE_DUNE_ISTL_GPU
        // GPU version based on dune-istl-gpu
        if( gpuSolverAdapter_ && assembledOperator_ )
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
            auto& matrix = assembledOperator_->matrixAdapter( *(gpuSolverAdapter_->parameter()) );

            // if preconditioner_ was set use that one, otherwise the one from the matrix object
            ISTLPreconditionerType& matrixPre = matrix.preconditionAdapter();
            ISTLPreconditionerType& precon    = ( preconditioner_ ) ? (*istlPre) : matrixPre;
            return solve( *gpuSolverAdapter_, matrix, scp, precon, u, w );
          }

          // only assembled implemented currently

          /*
          if constexpr (std::is_same< SolverDiscreteFunctionType, RangeFunctionType >::value )
          {
            assert( operator_ );
            assert( istlPre );
            typedef ISTLLinearOperatorAdapter< OperatorType >  ISTLOperatorType;
            ISTLOperatorType istlOperator( *operator_, space, space );
            return solve( *gpuSolverAdapter_, istlOperator, scp, *istlPre, u, w );
          }
          */
        }
        else
#endif
        {
          assert( solverAdapter_ );
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
            auto& matrix = assembledOperator_->matrixAdapter( *(solverAdapter_->parameter()) );
            // if preconditioner_ was set use that one, otherwise the one from the matrix object
            ISTLPreconditionerType& matrixPre = matrix.preconditionAdapter();
            ISTLPreconditionerType& precon    = ( preconditioner_ ) ? (*istlPre) : matrixPre;
            return solve( *solverAdapter_, matrix, scp, precon, u, w );
          }

          if constexpr (std::is_same< SolverDiscreteFunctionType, RangeFunctionType >::value )
          {
            assert( operator_ );
            assert( istlPre );
            typedef ISTLLinearOperatorAdapter< OperatorType >  ISTLOperatorType;
            ISTLOperatorType istlOperator( *operator_, space, space );
            return solve( *solverAdapter_, istlOperator, scp, *istlPre, u, w );
          }
        }

        DUNE_THROW(InvalidStateException,"ISTLInverseOperator::apply: No valid operator found!");
        return -1;
      }

      //! final solve execution only copying the right hand side
      template< class SolverAdapter, class OperatorAdapter, class ISTLPreconditioner, class DomainFunction >
      int solve ( SolverAdapter& solverAdapter, OperatorAdapter &istlOperator, ParallelScalarProductType &scp,
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
        solverAdapter.setMaxIterations( parameter_->maxIterations() );
        solverAdapter( istlOperator, scp, preconditioner, rhs_->blockVector(), w.blockVector(), result );
        return (result.converged) ? result.iterations : -(result.iterations);
      }

      using BaseType :: operator_;
      using BaseType :: assembledOperator_;
      using BaseType :: preconditioner_;

      using BaseType :: rhs_;
      using BaseType :: x_;

      using BaseType :: rightHandSideCopied_;
      using BaseType :: parameter_;

      mutable std::unique_ptr< SolverAdapterType > solverAdapter_;
      mutable std::unique_ptr< GPUSolverAdapterType > gpuSolverAdapter_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
