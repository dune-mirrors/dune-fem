#ifndef DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH
#define DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH

#include <limits>

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

      PreconditionedInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                      unsigned int maxIterations, bool verbose )
      : preconditioner_( op ),
        inverseOperator_( op, preconditioner_, redEps, absLimit, maxIterations, verbose )
      {}

      PreconditionedInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                      unsigned int maxIterations = std::numeric_limits< unsigned int >::max() )
      : preconditioner_( op ),
        inverseOperator_( op, preconditioner_, redEps, absLimit, maxIterations )
      {}

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        inverseOperator_( u, w );
      }

      unsigned int iterations () const { return inverseOperator_.iterations(); }

    private:
      Preconditioner preconditioner_;
      InverseOperator inverseOperator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_PRECONDITIONEDINVERSEOPERATOR_HH
