#ifndef DUNE_FEM_LOCALFUNCTION_INLINE_HH
#define DUNE_FEM_LOCALFUNCTION_INLINE_HH

#include "localfunction.hh"

namespace Dune
{

  // LocalFunctionDefault
  // --------------------

  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: operator+= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( 1, lf );
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: operator-= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( -1, lf );
  }

  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void 
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
  :: assign ( const LocalFunction< T > &lf ) 
  {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = lf[ i ];
  }

  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  inline void 
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
  :: clear ( ) 
  {
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = 0.0;
  }
  
  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const RangeFieldType s,
              const LocalFunction< T > &lf )
  {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] += s * lf[ i ];
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< int diffOrder, class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                  const PointType &x,
                  RangeType &ret ) const
  {
    asImp().baseFunctionSet().evaluate( diffVariable, x, asImp(), ret);
  }


  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: evaluate ( const PointType &x,
                  RangeType &ret ) const
  {
    asImp().baseFunctionSet().evaluate( x, asImp(), ret);
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: jacobian( const PointType &x,
                 JacobianRangeType &ret ) const
  {
    asImp().baseFunctionSet().jacobian( 
        x, 
        asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
        asImp(), ret);
  }
 

  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor )
  {
    asImp().baseFunctionSet().axpy( x, factor, asImp() );
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const JacobianRangeType &factor )
  {
    asImp().baseFunctionSet().axpy( 
            x, 
            asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
            factor, 
            asImp() );
  }

  
  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    asImp().baseFunctionSet().axpy( 
          x, 
          asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
          factor1, factor2, asImp() );
  }

} // end namespace Dune 
#endif
