namespace Dune
{

  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  void VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    :: evaluate ( const DomainType &x,
                  RangeType &ret ) const
  {
    const BaseFunctionSetType &baseFunctionSet = this->baseFunctionSet();

    ret = 0;

    const unsigned int numDofs = this->numDofs();
    for( unsigned int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseFunctionSet.evaluate( i, x, phi );
      ret.axpy( (*this)[ i ], phi );
    }
  }



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  template< class QuadratureType >
  void VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    :: evaluate ( const QuadratureType &quadrature,
                  const int quadPoint,
                  RangeType &ret ) const
  {
    const BaseFunctionSetType &baseFunctionSet = this->baseFunctionSet();

    ret = 0;

    const unsigned int numDofs = this->numDofs();
    for( unsigned int i = 0; i < numDofs; ++i )
    {
      RangeType phi;
      baseFunctionSet.evaluate( i, quadrature, quadPoint, phi );
      ret.axpy( (*this)[ i ], phi );
    }
  }




  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  void VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    :: init ( const Codim0EntityType &entity )
  {
    const DiscreteFunctionSpaceType &space = discreteFunction_.space();

    const bool multipleBaseSets
      = space.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      bool updateBaseSet = true;
      if( !multipleBaseSets && (geometry_ != 0) )
        updateBaseSet
          = (baseFunctionSet_.geometryType() != entity.geometry().type());

      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = space.baseFunctionSet( entity );

        const unsigned int numDofs = baseFunctionSet_.numBaseFunctions();
        values_.resize( numDofs );

        needCheckGeometry_ = space.multipleGeometryTypes();
      }
    }

    // cache geometry
    geometry_ = &(entity.geometry());

    assert( baseFunctionSet_.geometryType == geometry_->type() );
    const unsigned int size = values_.size();
    for( unsigned int i = 0; i < size; ++i )
    {
      const unsigned int index = space.mapToGlobal( entity, i );
      values_[ i ] = &(discreteFunction_.dofVector()[ index ]);
    }
  }



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  void VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    :: jacobian ( const DomainType &x,
                  JacobianRangeType &ret ) const
  {
    const BaseFunctionSetType &baseFunctionSet = this->baseFunctionSet();

    const GeometryJacobianInverseType &geoJacobianInverse
      = geometry_->jacobianInverseTransposed( x );

    ret = 0;
  
    const unsigned int numDofs = this->numDofs();
    for( unsigned int i = 0; i < numDofs; ++i )
    {
      JacobianRangeType DPhi;
      baseFunctionSet.jacobian( i, x, DPhi );
      ret.axpy( (*this)[ i ], DPhi );
    }

    for( unsigned int i = 0; i < DimRange; ++i )
      ret[ i ] = FMatrixHelp :: mult( geoJacobianInverse, ret[ i ] );
  }



  template< class DiscreteFunctionSpaceImp, class DofVectorImp >
  template< class QuadratureType >
  void VectorLocalFunction< DiscreteFunctionSpaceImp, DofVectorImp >
    :: jacobian ( const QuadratureType &quadrature,
                  const int quadPoint,
                  JacobianRangeType &ret ) const
  {
    const BaseFunctionSetType &baseFunctionSet = this->baseFunctionSet();

    const DomainType x = geometry_.global( quadrature.point[ quadPoint ] );

    const GeometryJacobianInverseType &geoJacobianInverse
      = geometry_->jacobianInverseTransposed( x );

    ret = 0;
  
    const unsigned int numDofs = this->numDofs();
    for( unsigned int i = 0; i < numDofs; ++i )
    {
      JacobianRangeType DPhi;
      baseFunctionSet.jacobian( i, quadrature, quadPoint, DPhi );
      ret.axpy( (*this)[ i ], DPhi );
    }

    for( unsigned int i = 0; i < DimRange; ++i )
      ret[ i ] = FMatrixHelp :: mult( geoJacobianInverse, ret[ i ] );
  }

}
