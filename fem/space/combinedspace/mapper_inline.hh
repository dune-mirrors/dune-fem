#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_INLINE_HH
#define DUNE_FEM_COMBINEDSPACE_MAPPER_INLINE_HH

#include "mapper.hh"

namespace Dune
{

  // CombinedMapper
  // --------------
  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline CombinedMapper< ContainedMapper, N, policy >
    :: CombinedMapper ( ContainedMapperType &mapper )
  : containedMapper_( &mapper ),
    utilLocal_( containedMapper(), numComponents ),
    utilGlobal_( containedMapper(), numComponents )
    //oldSize_( spc_.size() ),
    //size_( spc_.size() )
  {}


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy > :: size() const
  {
    return numComponents * containedMapper().size();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedMapper, N, policy >
    :: begin ( const EntityType &entity ) const
  {
    return DofMapIteratorType( containedMapper().begin( entity ), utilGlobal_ );
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedMapper, N, policy >
    :: end ( const EntityType &entity ) const
  {
    return DofMapIteratorType( containedMapper().end( entity ), utilGlobal_ );
  }

  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: mapToGlobal ( const EntityType &entity, const int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal
      = containedMapper().mapToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: mapEntityDofToGlobal ( const Entity &entity, int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal = 
      containedMapper().mapEntityDofToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: maxNumDofs () const
  {
    return numComponents * containedMapper().maxNumDofs();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numDofs ( const EntityType &entity ) const
  {
    return numComponents * containedMapper().numDofs( entity );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numEntityDofs ( const Entity &entity ) const
  {
    return numComponents * containedMapper().numEntityDofs( entity );
  }

  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: newSize () const
  {
    return numComponents * containedMapper().newSize();
  }
 

  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: newIndex ( const int hole, const int block ) const
  {
    if( policy == PointBased )
    {
      const int component = utilGlobal_.component( hole );
      const int containedHole = utilGlobal_.containedDof( hole );
      const int containedIndex = containedMapper().newIndex( containedHole, block );
      return utilGlobal_.combinedDof( containedIndex, component );
    }
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().newIndex( hole, containedBlock );
      return containedOffset + component * containedMapper().newSize();
    }
    
#if 0
    assert(policy != PointBased || block==0);
    return (policy == PointBased) ?
     (mapper_.newIndex(hole/N,block)*N+hole%N) :
     (mapper_.newIndex(hole,block%mapper_.numBlocks())+
      size_*block/mapper_.numBlocks());
    
    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.newSize(), Int2Type<policy>()));

    const int component = tmpUtilGlobal.component(hole);
    const int contained = tmpUtilGlobal.containedDof(hole);

    const int containedNew = mapper_.newIndex(contained,0);
    // const int containedNew = mapper_.newIndex(contained,block);

    return tmpUtilGlobal.combinedDof(containedNew, component);
#endif
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: oldIndex ( const int hole, const int block ) const
  {
    if( policy == PointBased )
    {
      const int component = utilGlobal_.component( hole );
      const int containedHole = utilGlobal_.containedDof( hole );
      const int containedIndex = containedMapper().oldIndex( containedHole, block );
      return utilGlobal_.combinedDof( containedIndex, component );
    }
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().oldIndex( hole, containedBlock );
      return containedOffset + component * containedMapper().size();
    }

#if 0
    assert(policy != PointBased || block==0);
    return (policy == PointBased) ?
     (mapper_.oldIndex(hole/N,0)*N+hole%N) :
     (mapper_.oldIndex(hole,block%mapper_.numBlocks())+
      oldSize_*block/mapper_.numBlocks());

    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.size(), Int2Type<policy>()));

    const int component = tmpUtilGlobal.component(hole);
    const int contained = tmpUtilGlobal.containedDof(hole);

    // const int containedNew = mapper_.oldIndex(contained,block);
    const int containedNew = mapper_.oldIndex(contained,0);

    return tmpUtilGlobal.combinedDof(containedNew, component);
#endif
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: numberOfHoles ( const int block ) const
  {
    if( policy == PointBased )
      return numComponents * containedMapper().numberOfHoles( block );
    else
      return containedMapper().numberOfHoles( block % containedMapper().numBlocks() );
  }


  
  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy > :: numBlocks () const
  {
    const int numContainedBlocks = containedMapper().numBlocks();
    if( policy == PointBased )
      return numContainedBlocks;
    else
      return numComponents * numContainedBlocks;
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline void CombinedMapper< ContainedMapper, N, policy > :: update ()
  {
    containedMapper().update();
    // calculate new size 
    //oldSize_ = size_;
    //size_ = containedMapper().size();
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: oldOffSet ( const int block ) const 
  {
    if( policy == PointBased )
      return numComponents * containedMapper().oldOffSet( block );
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().oldOffSet( containedBlock );
      return containedOffset + component * containedMapper().size();
    }
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedMapper, N, policy >
    :: offSet ( const int block ) const
  {
    if( policy == PointBased )
      return numComponents * containedMapper().offSet( block );
    else
    {
      const int numContainedBlocks = containedMapper().numBlocks();
      const int containedBlock = block % numContainedBlocks;
      const int component = block / numContainedBlocks;
      
      const int containedOffset = containedMapper().offSet( containedBlock );
      return containedOffset + component * containedMapper().newSize();
    }
  }



  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline bool CombinedMapper< ContainedMapper, N, policy >
    :: consecutive () const
  {
    return containedMapper().consecutive();
  }


  template< class ContainedMapper, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedMapper, N, policy > :: ContainedMapperType &
  CombinedMapper< ContainedMapper, N, policy > :: containedMapper () const
  {
    return *containedMapper_;
  }

  
#if 0
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: chooseSize ( int pointBased,
                    int variableBased,
                    Int2Type< PointBased > )
  {
    return pointBased;
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: chooseSize ( int pointBased,
                    int variableBased,
                    Int2Type< VariableBased > )
  {
    return variableBased;
  }
#endif

} // end namespace Dune

#endif
