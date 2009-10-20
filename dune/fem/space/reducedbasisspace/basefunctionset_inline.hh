#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH

#include "basefunctionset.hh"

namespace Dune
{

  template< class BaseFunctionImp >
  template< int diffOrd, class PointType >
  inline void ReducedBasisBaseFunctionSet< BaseFunctionImp >
    :: evaluate ( const int baseFunction,
                  const FieldVector< deriType, diffOrd > &diffVariable,
                  const PointType &x,
                  RangeType &phi ) const
  {
    assert( baseFunctionList_ != NULL );
    assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

    typedef typename LocalBaseFunctionType :: BaseFunctionSetType
      LocalBaseFunctionSetType;

    const LocalBaseFunctionType localBaseFunction
      = (*baseFunctionList_)[ baseFunction ]->localFunction( entity() );
    const LocalBaseFunctionSetType &localBaseFunctionSet
      = localBaseFunction.baseFunctionSet();
    const int numLocalBaseFunctions = localBaseFunctionSet.numBaseFunctions();
    
    phi = 0;
    for( int i = 0; i < numLocalBaseFunctions; ++i )
    {
      RangeType psi;
      localBaseFunctionSet.evaluate( i, diffVariable, x, psi );
      phi.axpy( localBaseFunction[ i ], psi );
    }
  }



  template< class BaseFunctionImp >
  template< class PointType >
  inline
  typename ReducedBasisBaseFunctionSet< BaseFunctionImp > :: RangeFieldType
  ReducedBasisBaseFunctionSet< BaseFunctionImp >
    :: evaluateSingle ( const int baseFunction,
                        const PointType &x,
                        RangeType &psi ) const
  {
    assert( baseFunctionList_ != NULL );
    assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

    typedef typename LocalBaseFunctionType :: BaseFunctionSetType
      LocalBaseFunctionSetType;

    const LocalBaseFunctionType localBaseFunction
      = (*baseFunctionList_)[ baseFunction ]->localFunction( entity() );
    const LocalBaseFunctionSetType &localBaseFunctionSet
      = localBaseFunction.baseFunctionSet();
    const int numLocalBaseFunctions = localBaseFunctionSet.numBaseFunctions();

    RangeFieldType ret = 0;
    for( int i = 0; i < numLocalBaseFunctions; ++i )
    {
      const RangeFieldType y
        = localBaseFunctionSet.evaluateSingle( i, x, psi );
      ret += localBaseFunction[ i ] * y;
    }
    return ret;
  }

}

#endif
