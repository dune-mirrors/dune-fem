namespace Dune 
{
  //- class SubSpace
  template <class CombinedSpaceImp>
  inline
  SubSpace<CombinedSpaceImp>::SubSpace(CombinedSpaceType& spc,
                                       int component) :
    BaseType(spc.gridPart()),
    spc_(spc),
    mapper_(spc, spc.mapper().containedMapper(), component),
    component_(component),
    baseSetMap_()
  {
    // create info for all geom types 
    AllGeomTypes<typename GridPartType :: IndexSetType, GridType>
      allGeomTypes(spc.gridPart().indexSet());

    // get types for codim 0 
    const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);

    //! this will not work for multiple base function sets 
    assert( spc.multipleBaseFunctionSets() == false );

    // create base sets for all existing geom types
    for(size_t i=0; i<geomTypes.size(); ++i)
    {
      if (baseSetMap_.find( geomTypes[i] ) == baseSetMap_.end()) 
      {
        baseSetMap_[ geomTypes[i] ]
          = new BaseFunctionSetImp( spc.baseFunctionSet( geomTypes[i] ) , component);
      }
    } // end for
  }
  
  template <class CombinedSpaceImp>
  inline
  SubSpace<CombinedSpaceImp>::~SubSpace()
  {
    typedef typename BaseFunctionMapType :: iterator iterator;
    iterator end = baseSetMap_.end();
    for (iterator it = baseSetMap_.begin(); it != end; ++it)
    {
      BaseFunctionSetImp * set = (BaseFunctionSetImp *) (*it).second;
      delete set; 
    }
  }
  
  //- class SubMapper
  template <class CombinedSpaceImp>
  inline
  int SubMapper<CombinedSpaceImp>::size() const 
  {
    return mapper_.size()/spc_.numComponents();
  }

  template <class CombinedSpaceImp>
  template <class EntityType>
  inline
  int SubMapper<CombinedSpaceImp>::
  mapToGlobal(EntityType& en, int localNum) const
  {
    utilGlobal_.newSize(mapper_.size()); // ok, since pointbased specialisation does nothing for newSize
    int globalNum = utilGlobal_.combinedDof(localNum, component_); 
    return mapper_.mapToGlobal(en,globalNum);
    
    const int containedGlobal = mapper_.mapToGlobal(en, localNum);

    return utilGlobal_.combinedDof(containedGlobal, component_);
  }

} // end namespace Dune
