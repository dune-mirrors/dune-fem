#ifndef DUNE_FEM_STANDARDLOCALFUNCTION_INLINE_HH
#define DUNE_FEM_STANDARDLOCALFUNCTION_INLINE_HH

#include "standardlocalfunction.hh"

namespace Dune
{

  // StandardLocalFunctionImpl
  // -------------------------

  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: StandardLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
  : discreteFunction_( discreteFunction ),
    values_( discreteFunction_.space().maxNumLocalDofs() ),
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
  inline int
  StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace > ::order () const
  {
    return discreteFunction_.space().order();
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
    :: init ( const EntityType &entity )
  {
    typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
    typedef typename BlockMapperType :: DofMapIteratorType DofMapIteratorType;
    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

    typedef typename DiscreteFunctionType :: DofBlockPtrType DofBlockPtrType;
    
    const DiscreteFunctionSpaceType &space = discreteFunction_.space();
    const bool multipleBaseSets = space.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      // if multiple base sets skip geometry call
      bool updateBaseSet = true;
      if( !multipleBaseSets && (entity_ != 0) )
        updateBaseSet = (baseFunctionSet_.geometryType() != entity.type());
      
      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = space.baseFunctionSet( entity );

        // note, do not use baseFunctionSet() here, entity might no have been set
        numDofs_ = baseFunctionSet_.numBaseFunctions();

        needCheckGeometry_ = space.multipleGeometryTypes();
      }
    }

    // cache entity
    entity_ = &entity;
    assert( baseFunctionSet_.geometryType() == entity.type() );

    assert( numDofs_ <= values_.size() );
    const BlockMapperType &mapper = space.blockMapper();
    const DofMapIteratorType end = mapper.end( entity );
    for( DofMapIteratorType it = mapper.begin( entity ); it != end; ++it )
    {
      assert( it.global() == mapper.mapToGlobal( entity, it.local() ) );
      
      DofBlockPtrType blockPtr = discreteFunction_.block( it.global() );
      
      const unsigned int localBlock = it.local() * blockSize;
      for( unsigned int i = 0; i < blockSize; ++i )
        values_[ localBlock + i ] = &((*blockPtr)[ i ]);
      
      // values_[ it.local() ] = &discreteFunction_.dof( it.global() );
    }
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int StandardLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: numDofs () const
  {
    return numDofs_;
  }

}
#endif
