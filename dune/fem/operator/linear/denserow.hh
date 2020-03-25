#ifndef DUNE_FEM_OPERATOR_LINEAR_DENSEROW_HH
#define DUNE_FEM_OPERATOR_LINEAR_DENSEROW_HH

#include <string>

#include <dune/common/dynmatrixev.hh>

#include <dune/fem/operator/matrix/densematrix.hh>

namespace Dune
{

  namespace Fem
  {

    // DenseRowLinearOperator
    // ----------------------

    template< class DomainFunction, class RangeFunction >
    class DenseRowLinearOperator
      : public DenseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
        public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef DenseRowLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef DenseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > BaseType;

    public:
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;

      static constexpr bool assembled = true;

      using BaseType::apply;

      DenseRowLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
        : BaseType( domainSpace, rangeSpace )
      {}

      virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const override { apply( arg, dest ); }

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_DENSEROW_HH
