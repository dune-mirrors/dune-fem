#ifndef DUNE_FEM_OPERATOR_LINEAR_ISTLADAPTER_HH
#define DUNE_FEM_OPERATOR_LINEAR_ISTLADAPTER_HH

#if HAVE_DUNE_ISTL

#include <dune/istl/operators.hh>

namespace Dune
{ 

  namespace Fem
  {

    // ISTLLinearOperatorAdapter
    // -------------------------

    template< class Operator >
    class ISTLLinearOperatorAdapter
      : public Dune::LinearOperator< typename Operator::DomainFunctionType::DofStorageType, typename Operator::RangeFunctionType::DofStorageType >
    {
      typedef ISTLLinearOperatorAdapter< Operator > ThisType;
      typedef Dune::LinearOperator< typename Operator::DomainFunctionType::DofStorageType, typename Operator::RangeFunctionType::DofStorageType > BaseType;

      typedef typename Operator::DomainFunctionType DomainFunctionType;
      typedef typename Operator::RangeFunctionType RangeFunctionType;

    public:
      enum {category=SolverCategory::sequential};

      typedef Operator OperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      typedef typename BaseType::domain_type domain_type;
      typedef typename BaseType::range_type range_type;
      typedef typename BaseType::field_type field_type;

      ISTLLinearOperatorAdapter ( const OperatorType &op,
                                  const DomainFunctionSpaceType &domainSpace,
                                  const RangeFunctionSpaceType &rangeSpace )
        : op_( op ),
          domainSpace_( domainSpace ),
          rangeSpace_( rangeSpace )
      {}

      virtual void apply ( const domain_type &x, range_type &y ) const
      {
        const DomainFunctionType fx( "ISTLLinearOperatorAdapter::apply::x", domainSpace_, x );
        RangeFunctionType fy( "ISTLLinearOperatorAdapter::apply::y", rangeSpace_, y );
        op_( fx, fy );
      }

      virtual void applyscaleadd ( field_type alpha, const domain_type &x, range_type &y ) const
      {
        const DomainFunctionType fx( "ISTLLinearOperatorAdapter::applyscaleadd::x", domainSpace_, x );
        RangeFunctionType fy( "ISTLLinearOperatorAdapter::applyscaleadd::y", rangeSpace_ );
        op_( fx, fy );
        y.axpy( alpha, fy.blockVector() );
      }

    private:
      const OperatorType &op_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };

  } // namespace Fem

} // namespace Dune 

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_ISTLADAPTER_HH
