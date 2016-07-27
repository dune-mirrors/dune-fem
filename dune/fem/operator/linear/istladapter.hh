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
      typedef typename Operator::RangeFunctionType  RangeFunctionType;

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

    protected:
      const OperatorType &op_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };

    template< class Operator >
    class ISTLMatrixFreeOperatorAdapter : public ISTLLinearOperatorAdapter< Operator >
    {
      typedef ISTLMatrixFreeOperatorAdapter< Operator >   ThisType;
      typedef ISTLLinearOperatorAdapter< Operator >       BaseType;

      typedef typename Operator::DomainFunctionType       DomainFunctionType;
      typedef typename Operator::RangeFunctionType        RangeFunctionType;

    public:
      enum {category=SolverCategory::sequential};

      typedef Operator OperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType  RangeFunctionSpaceType;

      static_assert( std::is_same< DomainFunctionType, ISTLBlockVectorDiscreteFunction< DomainFunctionSpaceType > > :: value,
                     "ISTLMatrixFreeOperatorAdapter only works with ISTLBlockVectorDiscreteFunction" );
      static_assert( std::is_same< RangeFunctionType,  ISTLBlockVectorDiscreteFunction< RangeFunctionSpaceType  > > :: value,
                     "ISTLMatrixFreeOperatorAdapter only works with ISTLBlockVectorDiscreteFunction" );

      typedef typename BaseType::domain_type domain_type;
      typedef typename BaseType::range_type range_type;
      typedef typename BaseType::field_type field_type;

      typedef Fem::ParallelScalarProduct< RangeFunctionType >                 ParallelScalarProductType;
      typedef Fem::IdentityPreconditionerWrapper< domain_type, range_type >   PreconditionAdapterType;

      ISTLMatrixFreeOperatorAdapter ( const OperatorType &op,
                                      const DomainFunctionSpaceType &domainSpace,
                                      const RangeFunctionSpaceType &rangeSpace )
        : BaseType( op, domainSpace, rangeSpace ),
          scp_( rangeSpace ),
          preconditioner_()
      {}

      //! return reference to preconditioner (here identity)
      PreconditionAdapterType& preconditionAdapter() { return preconditioner_; }
      //! return reference to preconditioner (here identity)
      const PreconditionAdapterType& preconditionAdapter() const { return preconditioner_; }

      virtual double residuum(const range_type& rhs, domain_type& x) const
      {
        range_type tmp( rhs );

        this->apply(x,tmp);
        tmp -= rhs;

        // return global sum of residuum
        return scp_.norm(tmp);
      }

      //! return reference to preconditioner
      ParallelScalarProductType& scp() const { return scp_; }

      //! dummy function returning 0
      double averageCommTime () const { return 0; }

    protected:
      mutable ParallelScalarProductType scp_;
      PreconditionAdapterType   preconditioner_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_ISTLADAPTER_HH
