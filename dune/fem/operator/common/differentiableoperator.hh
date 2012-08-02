#ifndef DUNE_FEM_DIFFERENTIABLEOPERATOR_HH
#define DUNE_FEM_DIFFERENTIABLEOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class DifferentiableOperator
     *  \brief abstract differentiable operator
     *
     *  Differentiable operators are operators providing a linearization.
     *
     *  \tparam  JacobianOperator  type of linear operator describing the Jacobian
     *                             (linearization) of this operator
     *
     *  \note The types for the operator's domain and range function are derived
     *        from the JacobianOperator.
     *
     *  \interfaceclass
     */
    template< class JacobianOperator >
    class DifferentiableOperator
    : public virtual Dune::Fem::Operator< typename JacobianOperator::DomainFunctionType,
                                          typename JacobianOperator::RangeFunctionType >
    {
      typedef Dune::Fem::Operator< typename JacobianOperator::DomainFunctionType,
                                   typename JacobianOperator::RangeFunctionType > BaseType;

    public:
      /** \brief type of linear operator modelling the operator's Jacobian */
      typedef JacobianOperator JacobianOperatorType;

      /** \brief type of discrete function in the operator's domain */
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      /** \brief type of discrete function in the operator's range */
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      /** \brief field type of the operator's domain */
      typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
      /** \brief field type of the operator's range */
      typedef typename RangeFunctionType::RangeFieldType RangeFieldType;

      /** \brief obtain linearization
       *
       *  \param[in]   u    argument discrete function
       *  \param[out]  jOp  destination Jacobian operator
       *
       *  \note This method has to be implemented by all derived classes.
       */
      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const = 0;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DIFFERENTIABLEOPERATOR_HH
