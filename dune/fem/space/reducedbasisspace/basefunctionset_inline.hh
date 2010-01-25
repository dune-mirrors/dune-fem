#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH

#include "basefunctionset.hh"

namespace Dune
{

  template< class BaseFunction >
  template< int diffOrd, class PointType >
  inline void ReducedBasisBaseFunctionSet< BaseFunction >
    ::evaluate ( const int baseFunction,
                 const FieldVector< deriType, diffOrd > &diffVariable,
                 const PointType &x,
                 RangeType &phi ) const
  {
    assert( baseFunctionList_ != 0 );
    assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

    const LocalBaseFunctionType localBaseFunction
      = (*baseFunctionList_)[ baseFunction ]->localFunction( entity() );
    localBaseFunction.evaluate( diffVariable, x, phi );
  }



  template< class BaseFunction >
  template< class PointType >
  inline typename ReducedBasisBaseFunctionSet< BaseFunction >::RangeFieldType
  ReducedBasisBaseFunctionSet< BaseFunction >
    ::evaluateSingle ( const int baseFunction,
                       const PointType &x,
                       RangeType &psi ) const
  {
    assert( baseFunctionList_ != 0 );
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

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH
