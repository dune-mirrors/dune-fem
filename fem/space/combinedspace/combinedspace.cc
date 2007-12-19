
namespace Dune {
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  CombinedSpace(GridPartType& gridpart) :
    BaseType(gridpart), 
    spc_(gridpart),
    mapper_(spc_, spc_.mapper()),
    baseSetMap_(),
    dm_(DofManagerFactoryType::getDofManager(spc_.grid()))
  {
    // get types for codim 0  
    AllGeomTypes<IndexSetType,typename Traits::GridType>
      allGeomTypes(gridpart.indexSet());
    
    const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);
    int maxNumDofs = -1;
    // create mappers and base sets for all existing geom types
    for(size_t i=0; i<geomTypes.size(); ++i)
    {
      if(baseSetMap_.find(geomTypes[i]) == baseSetMap_.end())
      {
        BaseFunctionSetImp* baseSet = 
          & SingletonProviderType::getObject(geomTypes[i]);
        // store in map 
        baseSetMap_[ geomTypes[i] ] = baseSet; 
        // calc max dofs 
        maxNumDofs = std::max(maxNumDofs,baseSet->numBaseFunctions());
      }
    }

  }
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  ~CombinedSpace() 
  {
    typedef typename BaseFunctionMapType :: iterator iterator;
    iterator end = baseSetMap_.end();
    for (iterator it = baseSetMap_.begin(); it != end; ++it)
    {
      BaseFunctionSetImp * set = (BaseFunctionSetImp *) (*it).second;
      SingletonProviderType::removeObject(*set);
    }
  }

  //- CombinedMapper
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::size() const 
  {
    return spc_.size()*N;
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class EntityType>
  inline int CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  mapToGlobal(EntityType& en, int localNum) const 
  {
    const int component = utilLocal_.component(localNum);
    const int containedLocal = utilLocal_.containedDof(localNum);
 
    const int containedGlobal = spc_.mapToGlobal(en, containedLocal);
    
    return utilGlobal_.combinedDof(containedGlobal, component);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  newIndex(const int hole, const int block) const 
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

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  oldIndex(const int hole, const int block) const 
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

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int  CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  numberOfHoles(const int block) const 
  {
    return (policy == PointBased) ?
     (mapper_.numberOfHoles(0)*N) : // in case of point based we have only on
     (mapper_.numberOfHoles(block%mapper_.numBlocks()));
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int  CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  numBlocks() const 
  {
    return (policy == PointBased) ? mapper_.numBlocks() : mapper_.numBlocks()*N;
  }

  //! update mapper, i.e. calculate new insertion points 
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline void CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  update()
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
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int  CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  oldOffSet(const int block) const 
  {
    return (policy == PointBased) ? 
      mapper_.oldOffSet(block)*N :
      mapper_.size()*(block/mapper_.numBlocks()) 
            + mapper_.oldOffSet(int(block/mapper_.numBlocks())); 
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int  CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  offSet(const int block) const 
  {
    return (policy == PointBased) ? 
      mapper_.offSet(block)*N :
      mapper_.newSize()*(block/mapper_.numBlocks()) 
            + mapper_.offSet(block/mapper_.numBlocks()); 
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline bool CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  needsCompress() const 
  {
    return mapper_.needsCompress ();
  }

} // end namespace Dune
