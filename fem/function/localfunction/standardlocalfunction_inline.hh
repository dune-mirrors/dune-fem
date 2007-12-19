namespace Dune
{

  // StandardLocalFunctionImpl
  // -------------------------

  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
  : discreteFunction_( discreteFunction ),
    values_(),
    baseFunctionSet_(),
    entity_( 0 ),
    numDofs_( 0 ),
    needCheckGeometry_( true )
  {
  }
 


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: StandardLocalFunctionImpl ( const ThisType &other )
  : discreteFunction_( other.discreteFunction_ ),
    values_( other.values_ ),
    baseFunctionSet_( other.baseFunctionSet_ ),
    entity_( other.entity_ ),
    numDofs_( other.numDofs_ ),
    needCheckGeometry_( other.needCheckGeometry_ )
  {
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: RangeFieldType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: operator[] ( const int num ) const
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }
  

  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: RangeFieldType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: operator[] ( const int num )
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: BaseFunctionSetType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: baseFunctionSet () const
  {
    assert( entity_ != 0 );
    return baseFunctionSet_;
  }


  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: EntityType &
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: entity () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline void StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
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



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: numDofs () const
  {
    return numDofs_;
  }



  // StandardLocalFunctionImpl (for CombinedSpace)
  // ---------------------------------------------

  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
  : discreteFunction_( discreteFunction ),
    values_(),
    baseFunctionSet_(),
    entity_( 0 ),
    numScalarDofs_( 0 ),
    needCheckGeometry_( true )
  {
  }



  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: StandardLocalFunctionImpl ( const ThisType &other )
  : discreteFunction_( other.discreteFunction_ ),
    values_( other.values_ ),
    baseFunctionSet_( other.baseFunctionSet_ ),
    entity_( other.entity_ ),
    numScalarDofs_( other.numScalarDofs_ ),
    needCheckGeometry_( other.needCheckGeometry_ )
  {
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline
  const typename StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: RangeFieldType &
  StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: operator[] ( const int num ) const
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline
  typename StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: RangeFieldType &
  StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: operator[] ( const int num )
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline
  const typename StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: BaseFunctionSetType &
  StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: baseFunctionSet () const
  {
    assert( entity_ != 0 );
    return baseFunctionSet_;
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline
  const typename StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: EntityType &
  StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: entity () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }


    
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline void StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
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


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline int StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: numDofs () const
  {
    return N * numScalarDofs_;
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline int StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: numScalarDofs () const
  {
    return numScalarDofs_;
  }


  
  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline void StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: mapLocalDofs ( const DofStoragePolicyType< PointBased > p,
                      const EntityType &entity,
                      const DiscreteFunctionSpaceType &space )
  {
    const int numDofs = this->numDofs();
    for( int i = 0; i < numDofs; i += N )
    {
      // const int baseindex = space.mapToGlobal( entity, i );
      int index = space.mapToGlobal( entity, i );
      for( int j = 0; j < N; ++j, ++index )
      {
        const int dof = i+j;
        // const int index = baseindex + j;
        assert( index == space.mapToGlobal( entity, dof ) );
        values_[ dof ] = &(discreteFunction_.dof( index ));
      }
    }
  }



  template< class DiscreteFunction,
            class ContainedFunctionSpace, int N, DofStoragePolicy policy >
  inline void StandardLocalFunctionImpl
    < DiscreteFunction, CombinedSpace< ContainedFunctionSpace, N, policy > >
    :: mapLocalDofs ( const DofStoragePolicyType< VariableBased > p,
                      const EntityType &entity,
                      const DiscreteFunctionSpaceType &space )
  {
    const int numDofs = this->numDofs();
    const int scalarSize = space.containedSpace().size();
    for( int i = 0; i < numDofs; i += N )
    {
      // const int baseindex = space.mapToGlobal( entity, i );
      int index = space.mapToGlobal( entity, i );
      for( int j = 0; j < N; ++j, index+=scalarSize )
      {
        const int dof = i + j;
        // const int index = baseindex + j * scalarSize;
        assert( index == space.mapToGlobal( entity, dof ) );
        values_[ dof ] = &(discreteFunction_.dof( index ));
      }
    }
  }
 
}
