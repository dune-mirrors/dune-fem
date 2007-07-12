#ifndef DUNE_FEM_IDENTITYOPERATOR_HH
#define DUNE_FEM_IDENTITYOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  template< class DiscreteFunctionImp >
  class IdentityOperator
  : public Operator< typename DiscreteFunctionImp :: RangeFieldType,
                     typename DiscreteFunctionImp :: RangeFieldType,
                     DiscreteFunctionImp,
                     DiscreteFunctionImp >
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;
    
    typedef DiscreteFunctionType DomainFunctionType;
    typedef DiscreteFunctionType RangeFunctionType;

    typedef typename DomainFunctionType :: RangeFieldType DomainFieldType;
    typedef typename RangeFunctionType :: RangeFieldType RangeFieldType;
   
  private:
    typedef IdentityOperator< DiscreteFunctionType > ThisType;
    typedef Operator< DomainFieldType, RangeFieldType, DomainFunctionType, RangeFunctionType >
      BaseType;

  public:
    inline IdentityOperator ()
    {
    }

    inline void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType;
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;

      ConstDofIteratorType uIt = u.dbegin();
      DofIteratorType wIt = w.dbegin();
      const DofIteratorType wEnd = w.dend();
      for( ; wIt != wEnd; ++wIt, ++uIt )
        *wIt = *uIt;
    }
  };
  
}

#endif
