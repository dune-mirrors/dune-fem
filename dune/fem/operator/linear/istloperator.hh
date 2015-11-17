#ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
#define DUNE_FEM_ISTLLINEAROPERATOR_HH

#if HAVE_DUNE_ISTL

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/istlmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! ISTLMatrixOperator
    template< class DomainFunction, class RangeFunction >
    struct ISTLLinearOperator
    : public ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef ISTLLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef ISTLMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;
      using BaseType::communicate;

      ISTLLinearOperator( const std::string & ,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace ) :
        BaseType( domainSpace, rangeSpace )
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
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
