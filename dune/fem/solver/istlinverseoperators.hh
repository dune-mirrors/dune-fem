#ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/operator/linear/istladapter.hh>

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
      typedef typename Preconditioner::RangeFunctionType RangeFunctionType;

    public:
      enum {category=SolverCategory::sequential};

      typedef typename BaseType::domain_type domain_type;
      typedef typename BaseType::range_type range_type;
      typedef typename BaseType::field_type field_type;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      ISTLPreconditionAdapter ( const Preconditioner *precon, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
      : precon_( precon ),
        domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace )
      {}

      // pre and post do nothing here
      virtual void pre ( domain_type &x, range_type &y ) {}
      virtual void post ( domain_type &x ) {}

      virtual void apply ( domain_type &x, const range_type &y )
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

    protected:
      const Preconditioner *precon_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };


    template< class BlockVector >
    struct ISTLSolverReduction
    {
      ISTLSolverReduction ( double redEps, double absLimit )
        : redEps_( redEps ),
          absLimit_( absLimit )
      {}

      double operator() ( const Dune::LinearOperator< BlockVector, BlockVector > &op,
                          Dune::ScalarProduct< BlockVector > &scp,
                          const BlockVector &rhs, const BlockVector &x ) const
      {
        if( absLimit_ < std::numeric_limits< double >::max() )
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
    };



    template< class Solver, class Reduction = ISTLSolverReduction< typename Solver::range_type > >
    struct ISTLSolverAdapter
    {
      typedef Solver SolverType;
      typedef Reduction ReductionType;

      typedef typename SolverType::domain_type domain_type;
      typedef typename SolverType::range_type range_type;

      ISTLSolverAdapter ( const ReductionType &reduction, unsigned int maxIterations, int verbose,
          const ParameterReader &parameter = Parameter::container() )
        : reduction_( reduction ),
          maxIterations_( maxIterations ),
          verbose_( verbose )
      {}

      template<class Op, class ScP, class PC >
      void operator () ( Op& op, ScP &scp, PC &pc,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
        int maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations_ );
        SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), maxIterations, verbose_ );
        solver.apply( x, rhs, result );
      }

    private:
      ReductionType reduction_;
      unsigned int maxIterations_;
      int verbose_;
    };


    template< class X, class Y, class F, class Reduction >
    struct ISTLSolverAdapter< Dune::RestartedGMResSolver< X, Y, F>, Reduction >
    {
      typedef Dune::RestartedGMResSolver< X, Y, F> SolverType;
      typedef Reduction ReductionType;

      typedef typename SolverType::domain_type domain_type;
      typedef typename SolverType::range_type range_type;

      ISTLSolverAdapter ( const ReductionType &reduction, unsigned int restart, unsigned int maxIterations, int verbose,
          const ParameterReader &parameter = Parameter::container() )
        : reduction_( reduction ),
          restart_( restart ),
          maxIterations_( maxIterations ),
          verbose_( verbose )
      {}

      ISTLSolverAdapter ( const Reduction &reduction, unsigned int maxIterations, int verbose,
          const ParameterReader &parameter = Parameter::container() )
        : reduction_( reduction ),
          restart_( parameter.getValue< int >( "fem.solver.gmres.restart", 20 ) ),
          maxIterations_( maxIterations ),
          verbose_( verbose )
      {}

      template<class Op, class ScP, class PC >
      void operator () ( Op& op, ScP &scp, PC &pc,
                         range_type &rhs, domain_type &x,
                         Dune::InverseOperatorResult &result ) const
      {
        int maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations_ );
        SolverType solver( op, scp, pc, reduction_( op, scp, rhs, x ), restart_, maxIterations, verbose_ );
        solver.apply( x, rhs, result );
      }

    private:
      ReductionType reduction_;
      unsigned int restart_;
      unsigned int maxIterations_;
      int verbose_;
    };


    template< class X >
    struct ISTLLoopSolver { typedef LoopSolver< X > Type; };

    template< class X >
    struct ISTLGradientSolver { typedef GradientSolver< X > Type; };

    template< class X >
    struct ISTLCGSolver { typedef CGSolver< X > Type; };

    template< class X >
    struct ISTLBiCGSTABSolver { typedef BiCGSTABSolver< X > Type; };

    template< class X >
    struct ISTLMINRESSolver { typedef MINRESSolver< X > Type; };

    template< class X >
    struct ISTLRestartedGMRes { typedef RestartedGMResSolver< X > Type; };


    // ISTLInverseOperator
    // -------------------

    template< class DiscreteFunction, template< class > class Solver,
              class Preconditioner = const Operator< DiscreteFunction, DiscreteFunction > >
    class ISTLInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Preconditioner PreconditionerType;

    protected:
      typedef ISTLLinearOperatorAdapter< OperatorType > ISTLOperatorType;
      typedef ISTLPreconditionAdapter< OperatorType > ISTLPreconditionerAdapterType;

      typedef Fem::ParallelScalarProduct< RangeFunctionType > ParallelScalarProductType;
      typedef typename DomainFunctionType::DofStorageType BlockVectorType;

      typedef ISTLSolverAdapter< typename Solver< BlockVectorType >::Type > SolverAdapterType;
      typedef typename SolverAdapterType::ReductionType ReductionType;
    public:

      typedef typename SolverAdapterType::SolverType SolverType;

      ISTLInverseOperator ( const OperatorType &op,
                            double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ISTLInverseOperator ( const OperatorType &op,
                            double redEps, double absLimit,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, nullptr, redEps, absLimit, std::numeric_limits< unsigned int >::max(), false, parameter ) {}

      ISTLInverseOperator ( const OperatorType &op,
                            double redEps, double absLimit, unsigned int maxIterations,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, nullptr, redEps, absLimit, maxIterations, false, parameter ) {}


      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(), false, parameter ) {}

      ISTLInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            double redEps, double absLimit, unsigned int maxIterations,
                            const ParameterReader &parameter = Parameter::container() )
        : ISTLInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, false, parameter ) {}


      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ISTLOperatorType istlOperator( operator_, w.space(), u.space() );
        ParallelScalarProductType scp( u.space() );

        if( !preconditioner_ )
        {
          ISTLPreconditionerAdapterType istlPreconditioner( nullptr, w.space(), u.space() );
          solve( istlOperator, scp, istlPreconditioner, u, w );
        }
        else
          solve( istlOperator, scp, *preconditioner_, u, w );
      }

      unsigned int iterations () const { return result_.iterations; }

    private:
      ISTLInverseOperator ( const OperatorType &op, PreconditionerType *preconditioner,
                            double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                            const ParameterReader &parameter )
        : operator_( op ),
          preconditioner_( preconditioner ),
          solverAdapter_( ReductionType( redEps, absLimit ), maxIterations, (Parameter::verbose() && verbose) ? 2 : 0, parameter )
      {}

      void solve ( ISTLOperatorType &istlOperator, ParallelScalarProductType &scp,
                   const OperatorType &preconditioner,
                   const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ISTLPreconditionerAdapterType istlPreconditioner( &preconditioner, w.space(), u.space() );
        solve( istlOperator, scp, istlPreconditioner, u, w );
      }

      template< class ISTLPreconditioner >
      void solve ( ISTLOperatorType &istlOperator, ParallelScalarProductType &scp,
                   ISTLPreconditioner &preconditioner,
                   const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        BlockVectorType rhs( u.blockVector() );
        solverAdapter_( istlOperator, scp, preconditioner, rhs, w.blockVector(), result_ );
      }

      const OperatorType &operator_;
      PreconditionerType *preconditioner_;
      SolverAdapterType solverAdapter_;
      mutable Dune::InverseOperatorResult result_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
