
namespace Dune {
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  CombinedSpace(GridPartType& gridpart) :
    BaseType(gridpart), 
    spc_(gridpart),
    mapper_(spc_, spc_.mapper()),
    baseSetVec_(GeometryIdentifier::numTypes, 0),
    subSpaces_(0),
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
      GeometryIdentifier::IdentifierType id =
                  GeometryIdentifier::fromGeo(geomTypes[i]);

      if(baseSetVec_[id] == 0 )
      {
        baseSetVec_[id] = & SingletonProviderType::getObject(geomTypes[i]);
        maxNumDofs = std::max(maxNumDofs,baseSetVec_[id]->numBaseFunctions());
      }
    }

    for (int i=0; i<N; ++i) 
      subSpaces_[i] = new SubSpaceType(*this,i);
  }
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  ~CombinedSpace() 
  {
    for (int i=0;i<N; ++i) delete subSpaces_[i];

    for (unsigned int i = 0; i < baseSetVec_.size(); ++i) 
    {
      if (baseSetVec_[i]) 
      {
        SingletonProviderType::removeObject(*baseSetVec_[i]);
        baseSetVec_[i] = 0;
      }
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
  newIndex(int num) const 
  {
    //assert( false );
    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.newSize(), Int2Type<policy>()));

    const int component = tmpUtilGlobal.component(num);
    const int contained = tmpUtilGlobal.containedDof(num);

    const int containedNew = mapper_.newIndex(contained);

    return tmpUtilGlobal.combinedDof(containedNew, component);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline int CombinedMapper<DiscreteFunctionSpaceImp, N, policy>::
  oldIndex(int num) const 
  {
    //assert(false);
    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.size(), Int2Type<policy>()));

    const int component = tmpUtilGlobal.component(num);
    const int contained = tmpUtilGlobal.containedDof(num);

    const int containedNew = mapper_.oldIndex(contained);

    return tmpUtilGlobal.combinedDof(containedNew, component);
  }
} // end namespace Dune
