#ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/solver/parameter.hh>

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
    :public Dune::Preconditioner< typename Preconditioner::RangeFunctionType::DofStorageType, typename Preconditioner::DomainFunctionType::DofStorageType >
    {
      typedef ISTLPreconditionAdapter< Preconditioner > ThisType;
      typedef Dune::Preconditioner< typename Preconditioner::RangeFunctionType::DofStorageType, typename Preconditioner::DomainFunctionType::DofStorageType > BaseType;

      typedef typename Preconditioner::DomainFunctionType DomainFunctionType;
      typedef typename Preconditioner::RangeFunctionType  RangeFunctionType;

    public:
#if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      enum {category=SolverCategory::sequential};
#endif // #if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)

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

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      SolverCategory::Category category () const override { return SolverCategory::sequential; }
#endif // #if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)

    protected:
      const Preconditioner *precon_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };


    template< class BlockVector >
    struct ISTLSolverReduction
    {
      ISTLSolverReduction ( double redEps, double absLimit, const SolverParameter& parameter )
        : redEps_( redEps ),
          absLimit_( absLimit ),
          errorMeasure_( parameter.errorMeasure() )
      {
      }

      double operator() ( const Dune::LinearOperator< BlockVector, BlockVector > &op,
                          Dune::ScalarProduct< BlockVector > &scp,
                          const BlockVector &rhs, const BlockVector &x ) const
      {

        if( errorMeasure_ == 0 && (absLimit_ < std::numeric_limits< double >::max()) )
        {
          BlockVector residuum( rhs );
          op.applyscaleadd( -1., x, residuum );
          const double res = scp.norm( residuum );
          return (res > 0 ? absLimit_ / res : 1e-3);
        }
        else
          return redEps_;
      }

    private:
      double redEps_;
      double absLimit_;
      int errorMeasure_ ;
    };

    template< int method,
              class X,
              class Reduction = ISTLSolverReduction< X > >
    struct ISTLSolverAdapter
    {
      typedef Reduction ReductionType;

      typedef X domain_type;
      typedef X range_type;

      ISTLSolverAdapter ( const ReductionType &reduction, unsigned int maxIterations, int verbose,
                          const SolverParameter& parameter )
        : reduction_( reduction ),
          method_( method < 0 ? parameter.krylovMethod() : method ),
          restart_( method_ == SolverParameter::gmres ? parameter.gmresRestart() : 0 ),
          maxIterations_( maxIterations ),
          verbose_( verbose ),
          parameter_( parameter )
      {
      }

      template<class Op, class ScP, class PC >
      void operator () ( Op& op, ScP &scp, PC &pc,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
        int maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations_ );
        if( method_ == SolverParameter::cg )
        {
          typedef Dune::CGSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::bicgstab )
        {
          typedef Dune::BiCGSTABSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::gmres )
        {
          typedef Dune::RestartedGMResSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), restart_, maxIterations, verbose_ );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::minres )
        {
          typedef Dune::MINRESSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
          solver.apply( x, rhs, result );
        }
        else if( method_ == SolverParameter::gradient )
        {
          typedef Dune::GradientSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
          solver.apply( x, rhs, result );
          return ;
        }
        else if( method_ == SolverParameter::loop )
        {
          typedef Dune::LoopSolver< X > SolverType;
          SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
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

      void setMaxIterations( unsigned int maxIterations ) { maxIterations_ = maxIterations; }

      const SolverParameter& parameter () const { return parameter_; }

    protected:
      template< class ImprovedMatrix >
      void callSuperLU ( ISTLParallelMatrixAdapterInterface< ImprovedMatrix >& op,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
#if HAVE_SUPERLU
        typedef typename ImprovedMatrix :: BaseType Matrix;
        const ImprovedMatrix& matrix = op.getmat();
        SuperLU< Matrix > solver( matrix, verbose_ );
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
      const unsigned int restart_;
      unsigned int maxIterations_;
      const int verbose_;

      SolverParameter parameter_;
    };



    // ISTLInverseOperator
    // -------------------

    template< class DiscreteFunction, int method = -1,
              class Preconditioner = const Operator< DiscreteFunction, DiscreteFunction > >
    class ISTLInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType  RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Preconditioner PreconditionerType;

    protected:
      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;

      typedef typename std::conditional<
         std::is_same< Dune::Fem::Operator< DomainFunctionType, RangeFunctionType >, OperatorType >::value,
                       Dune::Fem::ISTLLinearOperator< DomainFunctionType, RangeFunctionType >,
                       OperatorType > :: type  AssembledOperatorType;

      typedef ISTLLinearOperatorAdapter< OperatorType >       ISTLOperatorType;
      typedef ISTLPreconditionAdapter< PreconditionerType >   ISTLPreconditionerAdapterType;

      typedef typename RangeFunctionType :: ScalarProductType     ParallelScalarProductType;
      typedef typename RangeFunctionType::DofStorageType BlockVectorType;

      typedef ISTLSolverAdapter< method, BlockVectorType > SolverAdapterType;
      typedef typename SolverAdapterType::ReductionType ReductionType;
    public:
      ISTLInverseOperator ( double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader & parameter = Parameter::container() )
        : solverAdapter_( ReductionType( redEps, absLimit, SolverParameter(parameter) ), maxIterations, (Parameter::verbose() && verbose) ? 2 : 0, SolverParameter(parameter) )
      {}

      ISTLInverseOperator ( double redEps, double absLimit,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, std::numeric_limits< unsigned int >::max(), SolverParameter(parameter).verbose(), parameter ) {}

      ISTLInverseOperator ( double redEps, double absLimit,
                            const SolverParameter & parameter )
        : ISTLInverseOperator( redEps, absLimit, std::numeric_limits< unsigned int >::max(), parameter.verbose(), parameter.parameter() ) {}

      ISTLInverseOperator ( double redEps, double absLimit, unsigned int maxIterations,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, maxIterations, false, parameter ) {}

      ISTLInverseOperator ( const OperatorType &op,
                            double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator ( redEps, absLimit, maxIterations, verbose, parameter )
      {
        bind( op );
      }

      ISTLInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, parameter )
      {
        bind( op );
      }

      ISTLInverseOperator ( const OperatorType& op,
                            double redEps, double absLimit, unsigned int maxIterations,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, maxIterations, parameter )
      {
        bind( op );
      }

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, maxIterations, verbose, parameter )
      {
        bind( op, preconditioner );
      }

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, parameter )
      {
        bind( op, preconditioner );
      }

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit, unsigned int maxIterations,
                            const ParameterReader & parameter = Parameter::container() )
        : ISTLInverseOperator( redEps, absLimit, maxIterations, parameter )
      {
        bind( op, preconditioner );
      }

      void bind ( const OperatorType &op )
      {
        operator_ = &op;
        matrixOp_ = dynamic_cast<const AssembledOperatorType*>( &op );
      }

      void bind ( const OperatorType &op, PreconditionerType &preconditioner )
      {
        bind( op );
        preconditioner_ = &preconditioner;
      }

      void unbind () { operator_ = nullptr; matrixOp_ = nullptr; preconditioner_ = nullptr; rhs_.reset(); x_.reset(); }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        apply( u, w );
      }

      template <class DImpl, class RImpl>
      void operator() ( const DiscreteFunctionInterface< DImpl >&u,
                        DiscreteFunctionInterface< RImpl >& w ) const
      {
        apply( u, w );
      }

      unsigned int iterations () const { return result_.iterations; }
      void setMaxIterations ( unsigned int maxIterations ) { solverAdapter_.setMaxIterations( maxIterations ); }

    private:
      // apply for arbitrary domain function type and matching range function type
      template <class DomainFunction>
      void apply( const DomainFunction& u, RangeFunctionType& w ) const
      {
        auto& scp = w.scalarProduct();
        // u may not be a discrete function, therefore use w.space()
        const DiscreteFunctionSpaceType& space = w.space();
        ISTLPreconditionerAdapterType istlPreconditioner( preconditioner_, space, space );

        if( matrixOp_ )
        {
          ISTLMatrixParameter matparm( solverAdapter_.parameter().parameter() );
          auto& matrix = matrixOp_->matrixAdapter( matparm );
          // if preconditioner_ was set use that one, otherwise the one from the matrix object
          typedef Dune::Preconditioner< BlockVectorType, BlockVectorType > PreconditionerType;
          PreconditionerType& matrixPre = matrix.preconditionAdapter();
          PreconditionerType& istlPre   = istlPreconditioner;
          PreconditionerType& precon    = ( preconditioner_ ) ? istlPre : matrixPre;
          solve( matrix, scp, precon, u, w );
        }
        else
        {
          assert( operator_ );
          ISTLOperatorType istlOperator( *operator_, space, space );
          solve( istlOperator, scp, istlPreconditioner, u, w );
        }
      }

      // apply for arbitrary types of discrete function (only works if matrixOp_ is set)
      template <class DomainFunction, class RangeFunction>
      void apply( const DomainFunction& u, RangeFunction& w ) const
      {
        if( !matrixOp_ )
          DUNE_THROW(Dune::NotImplemented, "ISTLInverseOperator::operator() for matrix free operators only makes sense for fixed types of domain and range functions");

        if( ! x_ )
        {
          x_.reset( new RangeFunctionType( "ISTLInvOp::x", w.space() ) );
        }

        // copy right hand side since ISTL solvers seem to modify it
        x_->assign( w );

        apply( u, *x_ );

        // store result in destination
        w.assign( *x_ );
      }

      //! final solve execution only copying the right hand side
      template< class OperatorAdapter, class ISTLPreconditioner, class DomainFunction >
      void solve ( OperatorAdapter &istlOperator, ParallelScalarProductType &scp,
                   ISTLPreconditioner &preconditioner,
                   const DomainFunction& u,
                   RangeFunctionType& w ) const
      {
        if( ! rhs_ )
        {
          // u may not be a discrete function, therefore use w.space()
          rhs_.reset( new DomainFunctionType( "ISTLInvOp::rhs", w.space() ) );
        }

        // copy right hand side since ISTL solvers seem to modify it
        rhs_->assign( u );

        solverAdapter_( istlOperator, scp, preconditioner, rhs_->blockVector(), w.blockVector(), result_ );
      }

      mutable std::unique_ptr< DomainFunctionType > rhs_;
      mutable std::unique_ptr< RangeFunctionType  > x_;

      const OperatorType *operator_ = nullptr;
      const AssembledOperatorType* matrixOp_ = nullptr;
      PreconditionerType *preconditioner_ = nullptr;

      SolverAdapterType solverAdapter_;
      mutable Dune::InverseOperatorResult result_;
    };


    //////////////////////////////////////////////////////////////////////
    //  deprecated old types
    //////////////////////////////////////////////////////////////////////

    static const int ISTLLoopSolver     = SolverParameter :: loop ;
    static const int ISTLGradientSolver = SolverParameter :: gradient ;
    static const int ISTLCGSolver       = SolverParameter :: cg ;
    static const int ISTLBiCGSTABSolver = SolverParameter :: bicgstab ;
    static const int ISTLMINRESSolver   = SolverParameter :: minres ;
    static const int ISTLRestartedGMRes = SolverParameter :: gmres ;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLLoopOp = ISTLInverseOperator< DF, SolverParameter::loop >;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLMINResOp = ISTLInverseOperator< DF, SolverParameter::minres >;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLBICGSTABOp = ISTLInverseOperator< DF, SolverParameter::bicgstab >;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLGMResOp = ISTLInverseOperator< DF, SolverParameter::gmres >;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLCGOp = ISTLInverseOperator< DF, SolverParameter::cg >;

    template <class DF, class Op = Dune::Fem::Operator< DF, DF > >
    using ISTLSuperLU = ISTLInverseOperator< DF, SolverParameter::superlu >;

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
