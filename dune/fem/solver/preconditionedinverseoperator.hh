#ifndef DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH
#define DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH

#include <limits>
#include <memory>

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {

    // PreconditionedInverseOperator
    // -----------------------------

    template< class Preconditioner, class InverseOperator >
    class PreconditionedInverseOperator
    : public Operator< typename Preconditioner::RangeFunctionType >
    {
      typedef PreconditionedInverseOperator< Preconditioner, InverseOperator > ThisType;
      typedef Operator< typename Preconditioner::RangeFunctionType > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::DomainFunctionType RangeFunctionType;

      typedef typename InverseOperator::OperatorType OperatorType;

      PreconditionedInverseOperator ( double redEps, double absLimit, unsigned int maxIterations, bool verbose )
      : inverseOperator_( redEps, absLimit, maxIterations, verbose )
      {}

      PreconditionedInverseOperator ( double redEps, double absLimit,
                                      unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : inverseOperator_( redEps, absLimit, maxIterations )
      {}

      PreconditionedInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                      unsigned int maxIterations, bool verbose )
      : inverseOperator_( redEps, absLimit, maxIterations, verbose )
      {
        bind( op );
      }

      PreconditionedInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                      unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : inverseOperator_( redEps, absLimit, maxIterations )
      {
        bind( op );
      }

      void bind ( const OperatorType &op )
      {
        preconditioner_.reset( new Preconditioner( op ) );
        asssert( preconditioner_ );
        inverseOperator_.bind( op, *preconditioner_ );
      }
      void unbind() { inverseOperator_.unbind(); preconditioner_.reset(); }

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        assert( preconditioner_ );
        inverseOperator_( u, w );
      }

      unsigned int iterations () const { return inverseOperator_.iterations(); }
      void setMaxIterations ( unsigned int maxIterations ) const { inverseOperator_.setMaxIterations( maxIterations ); }

    private:
      InverseOperator inverseOperator_;
      std::unique_ptr< Preconditioner > preconditioner_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH
