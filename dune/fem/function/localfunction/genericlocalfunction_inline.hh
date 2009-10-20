#ifndef DUNE_FEM_GENERICLOCALFUNCTION_INLINE_HH
#define DUNE_FEM_GENERICLOCALFUNCTION_INLINE_HH

#include "genericlocalfunction.hh"

namespace Dune
{

  // GenericLocalFunctionImpl
  // -------------------------

  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::GenericLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
  : discreteFunction_( discreteFunction ),
    values_( discreteFunction_.space().maxNumLocalDofs() ),
    baseFunctionSet_(),
    entity_( 0 ),
    numDofs_( 0 ),
    needCheckGeometry_( true )
  {}
 


  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::GenericLocalFunctionImpl ( const ThisType &other )
  : discreteFunction_( other.discreteFunction_ ),
    values_( other.values_ ),
    baseFunctionSet_( other.baseFunctionSet_ ),
    entity_( other.entity_ ),
    numDofs_( other.numDofs_ ),
    needCheckGeometry_( other.needCheckGeometry_ )
  {}



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::RangeFieldType &
  GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::operator[] ( const int num ) const
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }
  

  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::RangeFieldType &
  GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::operator[] ( const int num )
  {
    assert( (num >= 0) && (num < numDofs()) );
    return *(values_[ num ]);
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int
  GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::order () const
  {
    return discreteFunction_.space().order();
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::BaseFunctionSetType &
  GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::baseFunctionSet () const
  {
    assert( entity_ != 0 );
    return baseFunctionSet_;
  }


  
  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline
  const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: EntityType &
  GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    :: entity () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline void GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::init ( const EntityType &entity )
  {
    typedef typename DiscreteFunctionSpaceType::MapperType MapperType;
   
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
    const MapperType &mapper = space.mapper();

    mapper.map( entity, indices_ );
    assert( indices_.size() == numDofs_ );
    for( unsigned int i = 0; i < numDofs_; ++i )
      values_[ i ] = &discreteFunction_.dof( indices_[ i ] );
  }



  template< class DiscreteFunction, class DiscreteFunctionSpace >
  inline int GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
    ::numDofs () const
  {
    return numDofs_;
  }

}

#endif // DUNE_FEM_GENERICLOCALFUNCTION_INLINE_HH
