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
  :: assign ( const LocalFunction< T > &lf ) {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = lf[ i ];
  }
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  inline void 
  LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
  :: clear ( ) {
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
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    ret = 0;
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, diffVariable, x, phi );
      ret.axpy( asImp()[ i ], phi );
    }
  }


  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: evaluate ( const PointType &x,
                  RangeType &ret ) const
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    ret = 0;
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      ret.axpy( asImp()[ i ], phi );
    }
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: jacobian( const PointType &x,
                 JacobianRangeType &ret ) const
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();
    
    JacobianRangeType refJacobian( 0 );
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        refJacobian[ j ].axpy( asImp()[ i ], grad[ j ] );
    }

    const GeometryJacobianInverseType &gjit
      = asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) );
    for( int i = 0; i < dimRange; ++i )
      // ret[ i ] = FMatrixHelp :: mult( gjit, refJacobian[ i ] );
      FieldMatrixHelper :: multiply( gjit, refJacobian[ i ], ret[ i ] );
  }
 

  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      asImp()[ i ] += phi * factor;
    }
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const JacobianRangeType &factor )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    JacobianRangeType factorInv;
    rightMultiply( factor, coordinate( x ), factorInv );
    
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        asImp()[ i ] += grad[ j ] * factorInv[ j ];
    }
  }

  
  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    JacobianRangeType factor2Inv;
    rightMultiply( factor2, coordinate( x ), factor2Inv );

    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      asImp()[ i ] += phi * factor1;
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        asImp()[ i ] += grad[ j ] * factor2Inv[ j ];
    }
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: rightMultiply( const JacobianRangeType &factor,
                      const DomainType &x,
                      JacobianRangeType &ret ) const
  {
    const GeometryJacobianInverseType &gjit
      = asImp().entity().geometry().jacobianInverseTransposed( x );

    FieldMatrixHelper :: multiply( factor, gjit, ret );
  }



  // LocalFunctionDefault (for CombinedSpace)
  // ----------------------------------------

  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: operator+= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( 1, lf );
  }



  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class T >
  inline void 
  LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: operator-= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( -1, lf );
  }


  
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class T >
  inline void 
  LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
  :: assign ( const LocalFunction< T > &lf ) {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = lf[ i ];
  }
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  inline void 
  LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
  :: clear ( ) {
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = 0.0;
  }
  
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: axpy ( const RangeFieldType s,
              const LocalFunction< T > &lf )
  {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] += s * lf[ i ];
  }


  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< int diffOrder, class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                  const PointType &x,
                  RangeType &ret ) const
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    ret = 0;
    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, diffVariable, x, phi );
      for( int j = 0; j < N; ++j )
        ret[ j ] += phi[ 0 ] * asImp()[ i*N + j ];
    }
  }


  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: evaluate ( const PointType &x,
                  RangeType &ret ) const
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    ret = 0;
    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      for( int j = 0; j < N; ++j )
        ret[ j ] += phi[ 0 ] * asImp()[ i*N + j ];
    }
  }



  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: jacobian( const PointType &x,
                 JacobianRangeType &ret ) const
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    const GeometryJacobianInverseType &gjit
      = asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) );

    ret = 0;

    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarJacobianRangeType gradPhiRef, gradPhi;
      baseSet.jacobianScalar( i, x, gradPhiRef );
      // gjit.umv( gradPhiRef[ 0 ], gradPhi[ 0 ] );
      FieldMatrixHelper :: multiply( gjit, gradPhiRef[ 0 ], gradPhi[ 0 ] );
      
      for( int j = 0; j < N; ++j )
        ret[ j ].axpy( asImp()[ i*N + j ], gradPhi[ 0 ] );
    }
  }


  
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      for( int j = 0; j < N; ++j )
        asImp()[ i*N + j ] += phi[ 0 ] * factor[ j ];
    }
  }

  
  
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const JacobianRangeType &factor )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();

    JacobianRangeType factorInv;
    rightMultiply( factor, coordinate( x ), factorInv );
    
    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, x, grad );
      for( int j = 0; j < N; ++j )
        asImp()[ i*N + j ] += grad[ 0 ] * factorInv[ j ];
    }
  }


  
  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: axpy ( const PointType &x,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    const BaseFunctionSetType &baseSet = asImp().baseFunctionSet();
    
    JacobianRangeType factor2Inv;
    rightMultiply( factor2, coordinate( x ), factor2Inv );
    
    const int numScalarDofs = asImp().numScalarDofs();
    for( int i = 0; i < numScalarDofs; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, x, grad );
      for( int j = 0; j < N; ++j )
        asImp()[ i*N + j ] += phi[ 0 ] * factor1[ j ] + grad[ 0 ] * factor2Inv[ j ];
    }
  }



  template< class ContainedFunctionSpace, int N, DofStoragePolicy policy,
            class LocalFunctionImp >
  inline void LocalFunctionDefault
    < CombinedSpace< ContainedFunctionSpace, N, policy >, LocalFunctionImp >
    :: rightMultiply( const JacobianRangeType &factor,
                      const DomainType &x,
                      JacobianRangeType &ret ) const
  {
    const GeometryJacobianInverseType &gjit
      = asImp().entity().geometry().jacobianInverseTransposed( x );

    FieldMatrixHelper :: multiply( factor, gjit, ret );
  }

}

#endif
