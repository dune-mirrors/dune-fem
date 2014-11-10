#ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
#define DUNE_FEM_ISTLLINEAROPERATOR_HH

#if HAVE_DUNE_ISTL

//- Dune fem includes 
#include <dune/fem/operator/matrix/istlmatrix.hh>

namespace Dune
{ 

  namespace Fem
  {

    // ISTLMatrixOperator
    // ------------------

    template< class DomainFunction, class RangeFunction >
    class ISTLLinearOperator
    : public ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef ISTLLinearOperator< DomainFunction, RangeFunction > This;
      typedef ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > Base;

    public:
      typedef typename Base::DomainSpaceType DomainSpaceType;
      typedef typename Base::RangeSpaceType RangeSpaceType;

      /** \copydoc Operator::assembled */
      static const bool assembled = true ;

      using Base::apply;
      using Base::communicate;

      ISTLLinearOperator ( const std::string &name,
                           const DomainSpaceType &domainSpace,
                           const RangeSpaceType &rangeSpace )
      : Base( domainSpace, rangeSpace )
      {}

      virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const
      {
        Base::apply( arg, dest );
      }

      const Base &systemMatrix () const
      {
        return *this;
      }
      Base &systemMatrix ()
      {
        return *this;
      }
    };

  } // namespace Fem

} // namespace Dune 

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
