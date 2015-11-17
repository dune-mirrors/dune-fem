#ifndef DUNE_FEM_SPOPERATOR_HH
#define DUNE_FEM_SPOPERATOR_HH

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! SparseRowLinearOperator
    template< class DomainFunction, class RangeFunction >
    struct SparseRowLinearOperator
    : public SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef SparseRowLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;

      SparseRowLinearOperator( const std::string & ,
                               const DomainSpaceType &domainSpace,
                               const RangeSpaceType &rangeSpace,
                               const std::string &paramfile = "" ) :
        BaseType( domainSpace, rangeSpace, paramfile )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      const BaseType &systemMatrix() const
      {
        return *this;
      }

      BaseType &systemMatrix()
      {
        return *this;
      }

      void communicate()
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPOPERATOR_HH
