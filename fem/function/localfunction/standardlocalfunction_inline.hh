namespace Dune
{

  // --------------------------------------------------------------------------
  // StandardLocalFunction
  // --------------------------------------------------------------------------

  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: init( const EntityType &entity )
  {
    const DiscreteFunctionSpaceType &space = discreteFunction_.space();
    const bool multipleBaseSets = space.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      // if multiple base sets skip geometry call
      bool updateBaseSet = true;
      if( !multipleBaseSets && (entity_ != 0) )
        updateBaseSet = (baseFunctionSet_.geometryType() != entity.geometry().type());
      
      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = space.baseFunctionSet( entity );

        // note, do not use baseFunctionSet() here, entity might no have been set
        numDofs_ = baseFunctionSet_.numBaseFunctions();
        values_.resize( numDofs_ );

        needCheckGeometry_ = space.multipleGeometryTypes();
      }
    }

    // cache entity
    entity_ = &entity;
    assert( baseFunctionSet_.geometryType() == entity.geometry().type() );

    for( int i = 0; i < numDofs_; ++i )
      values_[ i ] = &(discreteFunction_.dof( space.mapToGlobal( entity, i ) ));
  }


  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: evaluate ( const DomainType &x,
                  RangeType &ret ) const
  {
    ret = 0;
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      ret.axpy( (*values_[ i ]), phi );
    }
  }



  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  template< class QuadratureType >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: evaluate ( const QuadratureType &quadrature,
                  const int quadPoint,
                  RangeType &ret ) const
  {
    ret = 0;
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, quadrature, quadPoint, phi );
      ret.axpy( (*values_[ i ]), phi );
    }
  }


      
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: jacobian( const DomainType &x,
                 JacobianRangeType &ret ) const
  {
    JacobianRangeType refJacobian( 0 );
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        refJacobian[ j ].axpy( *(values_[ i ]), grad[ j ] );
    }

    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );
    for( int i = 0; i < dimRange; ++i )
      ret[ i ] = FMatrixHelp :: mult( gjit, refJacobian[ i ] );
  }


 
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  template< class QuadratureType >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: jacobian( const QuadratureType &quadrature,
                 const int quadPoint,
                 JacobianRangeType &ret ) const
  {
    const DomainType &x = quadrature.point( quadPoint );

    JacobianRangeType refJacobian( 0 );
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, quadrature, quadPoint, grad );
      for( int j = 0; j < dimRange; ++j )
        refJacobian[ j ].axpy( *(values_[ i ]), grad[ j ] );
    }

    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );
    for( int i = 0; i < dimRange; ++i )
      ret[ i ] = FMatrixHelp :: mult( gjit, refJacobian[ i ] );
  }


  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const DomainType &x,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      *(values_[ i ]) += phi * factor;
    }
  }

  
  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  template< class QuadratureType >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, quadrature, quadPoint, phi );
      *(values_[ i ]) += phi * factor;
    }
  }

  
  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const DomainType &x,
              const JacobianRangeType &factor )
  {
    JacobianRangeType factorInv;
    rightMultiply( factor, x, factorInv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        *(values_[ i ]) += grad[ j ] * factorInv[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  template< class QuadratureType >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const JacobianRangeType &factor )
  {
    JacobianRangeType factorInv;
    rightMultiply( factor, quadrature.point( quadPoint ), factorInv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      JacobianRangeType grad;
      baseSet.jacobian( i, quadrature, quadPoint, grad );
      for( int j = 0; j < dimRange; ++j )
        *(values_[ i ]) += grad[ j ] * factorInv[ j ];
    }
  }


  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const DomainType &x,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    JacobianRangeType factor2Inv;
    rightMultiply( factor2, x, factor2Inv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, x, phi );
      *(values_[ i ]) += phi * factor1;
      JacobianRangeType grad;
      baseSet.jacobian( i, x, grad );
      for( int j = 0; j < dimRange; ++j )
        *(values_[ i ]) += grad[ j ] * factor2Inv[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  template< class QuadratureType >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    JacobianRangeType factor2Inv;
    rightMultiply( factor2, quadrature.point( quadPoint ), factor2Inv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numDofs_; ++i )
    {
      RangeType phi;
      baseSet.evaluate( i, quadrature, quadPoint, phi );
      *(values_[ i ]) += phi * factor1;
      JacobianRangeType grad;
      baseSet.jacobian( i, quadrature, quadPoint, grad );
      for( int j = 0; j < dimRange; ++j )
        *(values_[ i ]) += grad[ j ] * factor2Inv[ j ];
    }
  }


    
  template< class DiscreteFunctionImp, class DiscreteFunctionSpaceImp >
  inline void StandardLocalFunction< DiscreteFunctionImp, DiscreteFunctionSpaceImp >
    :: rightMultiply( const JacobianRangeType &factor,
                      const DomainType &x,
                      JacobianRangeType &ret ) 
  {
    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );
    
    const int rows = JacobianRangeType :: rows;
    for( int i = 0; i < rows; ++i )
    {
      const int cols = JacobianRangeType :: cols;
      for( int j = 0; j < cols; ++j )
      {
        RangeFieldType value( 0 );
        for( int k = 0; k < cols; ++k )
          value += factor[ i ][ k ] * gjit[ k ][ j ];
        ret[ i ][ j ] = value;
      }
    }
  }



  // --------------------------------------------------------------------------
  // StandardLocalFunction (for CombinedSpace)
  // --------------------------------------------------------------------------
 
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: init( const EntityType &entity )
  {
    const DiscreteFunctionSpaceType &space = discreteFunction_.space();
    const bool multipleBaseSets = space.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      // if multiple base sets skip geometry call
      bool updateBaseSet = true;
      if( !multipleBaseSets && (entity_ != 0) )
        updateBaseSet = (baseFunctionSet_.geometryType() != entity.geometry().type());
      
      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = space.baseFunctionSet( entity );

        // note, do not use baseFunctionSet() here, entity might no have been set
        numScalarDofs_ = baseFunctionSet_.numDifferentBaseFunctions();
        values_.resize( N * numScalarDofs_ );

        needCheckGeometry_ = space.multipleGeometryTypes();
      }
    }

    // cache entity
    entity_ = &entity;
    assert( baseFunctionSet_.geometryType() == entity.geometry().type() );

    assert( values_.size() == (unsigned int)(N * numScalarDofs_) );
    DofStoragePolicyType< policy > p;
    mapLocalDofs( p, entity, space );
  }



  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: evaluate ( const DomainType &x,
                  RangeType &ret ) const
  {
    ret = 0;
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      for( int j = 0; j < N; ++j )
        ret[ j ] += phi[ 0 ] * (*(values_[ i*N + j ]));
    }
  }


  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  template< class QuadratureType >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: evaluate ( const QuadratureType &quadrature,
                  const int quadPoint,
                  RangeType &ret ) const
  {
    ret = 0;
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, quadrature, quadPoint, phi );
      for( int j = 0; j < N; ++j )
        ret[ j ] += phi[ 0 ] * (*(values_[ i*N + j ]));
    }
  }


      
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: jacobian( const DomainType &x,
                 JacobianRangeType &ret ) const
  {
    ret = 0;

    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );

    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarJacobianRangeType gradPhiRef, gradPhi;
      baseSet.jacobianScalar( i, x, gradPhiRef );
      gjit.umv( gradPhiRef[ 0 ], gradPhi[ 0 ] );
      
      for( int j = 0; j < N; ++j )
        ret[ j ].axpy( *(values_[ i*N + j ]), gradPhi[ 0 ] );
    }
  }


 
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  template< class QuadratureType >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: jacobian( const QuadratureType &quadrature,
                 const int quadPoint,
                 JacobianRangeType &ret ) const
  {
    ret = 0;
     
    const DomainType &x = quadrature.point( quadPoint );
    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );

    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarJacobianRangeType gradPhiRef, gradPhi;
      baseSet.jacobianScalar( i, quadrature, quadPoint, gradPhiRef );
      gjit.umv( gradPhiRef[ 0 ], gradPhi[ 0 ] );
      
      for( int j = 0; j < N; ++j )
        ret[ j ].axpy( *(values_[ i*N + j ]), gradPhi[ 0 ] );
    }
  }


  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const DomainType &x,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += phi[ 0 ] * factor[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  template< class QuadratureType >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const RangeType &factor )
  {
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, quadrature, quadPoint, phi );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += phi[ 0 ] * factor[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const DomainType &x,
              const JacobianRangeType &factor )
  {
    JacobianRangeType factorInv;
    rightMultiply( factor, x, factorInv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, x, grad );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += grad[ 0 ] * factorInv[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  template< class QuadratureType >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const JacobianRangeType &factor )
  {
    JacobianRangeType factorInv;
    rightMultiply( factor, quadrature.point( quadPoint ), factorInv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, quadrature, quadPoint, grad );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += grad[ 0 ] * factorInv[ j ];
    }
  }


  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const DomainType &x,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    JacobianRangeType factor2Inv;
    rightMultiply( factor2, x, factor2Inv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, x, phi );
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, x, grad );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += phi[ 0 ] * factor1[ j ] + grad[ 0 ] * factor2Inv[ j ];
    }
  }

  
  
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  template< class QuadratureType >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: axpy ( const QuadratureType &quadrature,
              const int quadPoint,
              const RangeType &factor1,
              const JacobianRangeType &factor2 )
  {
    JacobianRangeType factor2Inv;
    rightMultiply( factor2, quadrature.point( quadPoint ), factor2Inv );
    
    const BaseFunctionSetType &baseSet = baseFunctionSet();
    for( int i = 0; i < numScalarDofs_; ++i )
    {
      ScalarRangeType phi;
      baseSet.evaluateScalar( i, quadrature, quadPoint, phi );
      ScalarJacobianRangeType grad;
      baseSet.jacobianScalar( i, quadrature, quadPoint, grad );
      for( int j = 0; j < N; ++j )
        *(values_[ i*N + j ]) += phi[ 0 ] * factor1[ j ] + grad[ 0 ] * factor2Inv[ j ];
    }
  }


    
  template< class DiscreteFunctionImp,
            class ContainedFunctionSpaceImp, int N, DofStoragePolicy policy >
  inline void StandardLocalFunction
    < DiscreteFunctionImp, CombinedSpace< ContainedFunctionSpaceImp, N, policy > >
    :: rightMultiply( const JacobianRangeType &factor,
                      const DomainType &x,
                      JacobianRangeType &ret )
  {
    const GeometryJacobianInverseType &gjit
      = entity().geometry().jacobianInverseTransposed( x );
    
    const int rows = JacobianRangeType :: rows;
    for( int i = 0; i < rows; ++i )
    {
      const int cols = JacobianRangeType :: cols;
      for( int j = 0; j < cols; ++j )
      {
        RangeFieldType value( 0 );
        for( int k = 0; k < cols; ++k )
          value += factor[ i ][ k ] * gjit[ k ][ j ];
        ret[ i ][ j ] = value;
      }
    }
  }
  
}
