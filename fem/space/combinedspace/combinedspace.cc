namespace Dune
{
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>
    :: CombinedSpace( GridPartType &gridpart )
  : BaseType( gridpart ),
    containedSpace_( gridpart ),
    mapper_( containedSpace_.mapper() ),
    blockMapper_( Traits :: BlockTraits :: containedBlockMapper( containedSpace_ ) ),
    baseSetMap_(),
    dm_( DofManagerFactoryType :: getDofManager( containedSpace_.grid() ) )
  {
    const std::vector<GeometryType>& geomTypes = containedSpace_.geomTypes(0);
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

} // end namespace Dune
