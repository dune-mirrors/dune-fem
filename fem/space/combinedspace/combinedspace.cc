namespace Dune
{
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline CombinedSpace<DiscreteFunctionSpaceImp, N, policy>
    :: CombinedSpace( GridPartType &gridpart )
  : BaseType( gridpart ),
    spc_( gridpart ),
    mapper_( spc_.mapper() ),
    blockMapper_( spc_.blockMapper() ),
    baseSetMap_(),
    dm_( DofManagerFactoryType :: getDofManager( spc_.grid() ) )
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

} // end namespace Dune
