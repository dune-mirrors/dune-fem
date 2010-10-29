#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH

#include "basefunctionset.hh"

namespace Dune
{

  template< class BaseFunction >
  template< int diffOrd, class PointType >
  inline void ReducedBasisBaseFunctionSet< BaseFunction >
    ::evaluate ( const int baseFunction,
                 const FieldVector< int, diffOrd > &diffVariable,
                 const PointType &x,
                 RangeType &phi ) const
  {
    assert( baseFunctionList_ != 0 );
    assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

    const LocalBaseFunctionType localBaseFunction
      = (*baseFunctionList_)[ baseFunction ]->localFunction( entity() );
    localBaseFunction.evaluate( diffVariable, x, phi );
  }
}

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_INLINE_HH
