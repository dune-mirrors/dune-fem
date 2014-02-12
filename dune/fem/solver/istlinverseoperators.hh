#ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH

#include <dune/common/nullptr.hh>

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/operator/linear/istladapter.hh>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioner.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declartion
    template< class DiscreteFunction >
    class ISTLGeneralizedMinResInverseOperator;

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
      typedef typename BaseType::domain_type domain_type;
      typedef typename BaseType::range_type range_type;
      typedef typename BaseType::field_type field_type;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      ISTLPreconditionAdapter ( Preconditioner *precon, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
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
      Preconditioner *precon_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };


    // ISLTGeneralizedMinResInverseOperator
    // ------------------------------------

    template< class DiscreteFunction >
    class ISTLGeneralizedMinResInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;


    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

    protected:
      typedef ISTLLinearOperatorAdapter< OperatorType > ISTLOperatorType;
      typedef ISTLPreconditionAdapter< PreconditionerType > ISTLPreconditionerType;

      typedef ParallelScalarProduct< typename RangeFunctionType::DiscreteFunctionSpaceType > ParallelScalarProductType;

      typedef typename DomainFunctionType::DofStorageType BlockVectorType;
    public:

      ISTLGeneralizedMinResInverseOperator ( const OperatorType &op,
                                             double redEps, double absLimit,
                                             unsigned int maxIterations, bool verbose )
      : operator_( op ),
        preconditioner_( nullptr ),
        redEps_( redEps ),
        absLimit_( absLimit ),
        maxIterations_( maxIterations ),
        verbose_( verbose )
      {
      }

      ISTLGeneralizedMinResInverseOperator ( const OperatorType &op,
                                             double redEps, double absLimit,
                                             unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : operator_( op ),
        preconditioner_( nullptr ),
        redEps_( redEps ),
        absLimit_( absLimit ),
        maxIterations_( maxIterations ),
        verbose_( false )
      {
      }

      ISTLGeneralizedMinResInverseOperator ( const OperatorType &op,
                                             const PreconditionerType &preconditioner,
                                             double redEps, double absLimit,
                                             unsigned int maxIterations, bool verbose )
      : operator_( op ),
        preconditioner_( &preconditioner ),
        redEps_( redEps ),
        absLimit_( absLimit ),
        maxIterations_( maxIterations ),
        verbose_( verbose )
      {
      }

      ISTLGeneralizedMinResInverseOperator ( const OperatorType &op,
                                             const PreconditionerType &preconditioner,
                                             double redEps, double absLimit,
                                             unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : operator_( op ),
        preconditioner_( &preconditioner ),
        redEps_( redEps ),
        absLimit_( absLimit ),
        maxIterations_( maxIterations ),
        verbose_( false )
      {
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ISTLOperatorType istlOperator( operator_, w.space(), u.space() );
        ISTLPreconditionerType istlPreconditioner( preconditioner_, w.space(), u.space() );
        ParallelScalarProductType scp( u.space() );

        // verbose only in verbose mode and for rank 0 
        int verb = ( Parameter :: verbose() && verbose_ ) ? 2 : 0;

        double reduction = redEps_;
        if( absLimit_ < std::numeric_limits< double >::max() )
        {
          BlockVectorType residuum( u.blockVector() );
          istlOperator.applyscaleadd( -1., w.blockVector(), residuum );
          const double res = scp.norm( residuum );
          reduction = (res > 0) ? absLimit_/ res : 1e-3;
        }

        RestartedGMResSolver< BlockVectorType > solver( istlOperator, scp, istlPreconditioner, reduction, paramRestart(), maxIterations_, verb );

        BlockVectorType rhs( u.blockVector() );
        InverseOperatorResult returnInfo;
        solver( rhs, w.blockVector(), returnInfo );

        iterations_ = returnInfo.iterations;
      }

      unsigned int iterations () const
      {
        return iterations_;
      }

    private:
      static int paramRestart ()
      {
        return Parameter::getValue< int >( "fem.solver.gmres.restart", 20 );
      }

      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
      const double redEps_;
      const double absLimit_;
      const int maxIterations_;
      const bool verbose_;
      mutable unsigned int iterations_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTLINVERSEOPERATORS_HH
