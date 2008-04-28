#ifndef DUNE_FEM_TEMPORARYLOCALFUNCTION_INLINE_HH
#define DUNE_FEM_TEMPORARYLOCALFUNCTION_INLINE_HH

#include "temporarylocalfunction.hh"

namespace Dune
{

  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: TemporaryLocalFunctionImpl ( const DiscreteFunctionSpaceType &dfSpace )
  : discreteFunctionSpace_( dfSpace ),
    entity_( 0 ),
    baseFunctionSet_(),
    dofs_( discreteFunctionSpace_.maxNumLocalDofs() ),
    needCheckGeometry_( true )
  {
  }


  
  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: TemporaryLocalFunctionImpl ( const DiscreteFunctionSpaceType &dfSpace,
                                    const EntityType &entity )
  : discreteFunctionSpace_( dfSpace ),
    entity_( &entity ),
    baseFunctionSet_( discreteFunctionSpace_.baseFunctionSet( entity ) ),
    dofs_( discreteFunctionSpace_.maxNumLocalDofs() ),
    needCheckGeometry_( true )
  {
  }



  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: TemporaryLocalFunctionImpl ( const ThisType &other )
  : discreteFunctionSpace_( other.discreteFunctionSpace_ ),
    entity_( other.entity_ ),
    baseFunctionSet_( other.baseFunctionSet_ ),
    dofs_( other.dofs_ ),
    needCheckGeometry_( true )
  {
  }


  
  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline
  const typename TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: RangeFieldType &
  TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: operator[] ( const int num ) const
  {
    assert( num < numDofs() );
    return dofs_[ num ];
  }



  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline
  typename TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: RangeFieldType &
  TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: operator[] ( const int num )
  {
    assert( num < numDofs() );
    return dofs_[ num ];
  }



  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline
  const typename TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: BaseFunctionSetType &
  TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: baseFunctionSet () const
  {
    assert( entity_ != 0 );
    return baseFunctionSet_;
  }


  
  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline
  const typename TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: EntityType &
  TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: entity () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }



  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline void TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: init ( const EntityType &entity )
  {
    const bool multipleBaseSets = discreteFunctionSpace_.multipleBaseFunctionSets();

    if( multipleBaseSets || needCheckGeometry_ )
    {
      // if multiple base sets skip geometry call
      bool updateBaseSet = true;
      if( !multipleBaseSets && (entity_ != 0) )
        updateBaseSet = (baseFunctionSet_.geometryType() != entity.geometry().type());
      
      if( multipleBaseSets || updateBaseSet )
      {
        baseFunctionSet_ = discreteFunctionSpace_.baseFunctionSet( entity );
        needCheckGeometry_ = discreteFunctionSpace_.multipleGeometryTypes();
      }
    }
    assert( (unsigned int)baseFunctionSet_.numBaseFunctions() <= dofs_.size() );

    entity_ = &entity;
    assert( baseFunctionSet_.geometryType() == entity.geometry().type() );
  }



  template< class DiscreteFunctionSpace, template< class > class ArrayAllocator >
  inline int TemporaryLocalFunctionImpl< DiscreteFunctionSpace, ArrayAllocator >
    :: numDofs () const
  {
    return baseFunctionSet_.numBaseFunctions();
  }

}

#endif
