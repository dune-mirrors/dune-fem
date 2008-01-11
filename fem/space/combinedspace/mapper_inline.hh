namespace Dune
{

  // CombinedMapper
  // --------------
  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline CombinedMapper< ContainedSpace, N, policy >
    :: CombinedMapper ( const ContainedDiscreteFunctionSpaceType &spc,
                        ContainedMapperType &mapper )
  : spc_( spc ),
    mapper_( mapper ),
    utilLocal_( spc_, numComponents ),
    utilGlobal_( spc_, numComponents ),
    oldSize_( spc_.size() ),
    size_( spc_.size() )
  {}



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy > :: size() const
  {
    return spc_.size() * numComponents;
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedSpace, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedSpace, N, policy >
    :: begin ( const EntityType &entity ) const
  {
    return DofMapIteratorType( mapper_.begin( entity ), utilGlobal_ );
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedSpace, N, policy > ::  DofMapIteratorType
  CombinedMapper< ContainedSpace, N, policy >
    :: end ( const EntityType &entity ) const
  {
    return DofMapIteratorType( mapper_.end( entity ), utilGlobal_ );
  }

  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: mapToGlobal ( const EntityType &entity, int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal = spc_.mapToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }

  template< class ContainedSpace, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: mapEntityDofToGlobal ( const Entity &entity, int localDof ) const 
  {
    const int component = utilLocal_.component( localDof );
    const int containedLocal = utilLocal_.containedDof( localDof );
 
    const int containedGlobal = 
      spc_.mapper().mapEntityDofToGlobal( entity, containedLocal );
    
    return utilGlobal_.combinedDof( containedGlobal, component );
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: maxNumDofs () const
  {
    return mapper_.maxNumDofs() * numComponents;
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: numDofs ( const EntityType &entity ) const
  {
    return mapper_.numDofs( entity ) * numComponents;
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  template< class Entity >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: numEntityDofs ( const Entity &entity ) const
  {
    return mapper_.numEntityDofs( entity ) * numComponents;
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: newSize () const
  {
    return mapper_.newSize() * numComponents;
  }
 

  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: newIndex ( const int hole, const int block ) const
  {
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
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: oldIndex ( const int hole, const int block ) const
  {
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
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: numberOfHoles ( const int block ) const
  {
    return (policy == PointBased) ?
     (mapper_.numberOfHoles(0)*N) : // in case of point based we have only on
     (mapper_.numberOfHoles(block%mapper_.numBlocks()));
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy > :: numBlocks () const
  {
    return (policy == PointBased) ? mapper_.numBlocks() : mapper_.numBlocks()*N;
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline void CombinedMapper< ContainedSpace, N, policy > :: update ()
  {
    // assert(0);
    // assure that update is only called once per 
    // dof manager resize 
    if( 1 ) // sequence_ != dm_.sequence() )
    {
      // calculate new size 
      oldSize_ = size_;
      size_ = spc_.size();
      // sequence_ = dm_.sequence();
    }
    else 
    {
      DUNE_THROW(InvalidStateException,"update of mapper should only be called once per sequence");
    }
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: oldOffSet ( const int block ) const 
  {
    return (policy == PointBased) ? 
      mapper_.oldOffSet(block)*N :
      mapper_.size()*(block/mapper_.numBlocks()) 
            + mapper_.oldOffSet(int(block/mapper_.numBlocks())); 
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline int CombinedMapper< ContainedSpace, N, policy >
    :: offSet ( const int block ) const
  {
    return (policy == PointBased) ? 
      mapper_.offSet(block)*N :
      mapper_.newSize()*(block/mapper_.numBlocks()) 
            + mapper_.offSet(block/mapper_.numBlocks()); 
  }



  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline bool CombinedMapper< ContainedSpace, N, policy >
    :: needsCompress () const
  {
    return mapper_.needsCompress ();
  }


  
  template< class ContainedSpace, int N, DofStoragePolicy policy >
  inline typename CombinedMapper< ContainedSpace, N, policy > :: ContainedMapperType &
  CombinedMapper< ContainedSpace, N, policy > :: containedMapper () const
  {
    return mapper_;
  }


  
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

} // end namespace Dune
