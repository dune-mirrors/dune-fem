namespace Dune 
{
  //- class SubSpace
  template <class CombinedSpaceImp>
  inline
  SubSpace<CombinedSpaceImp>::SubSpace(const CombinedSpaceType& spc,
                                       int component) :
    BaseType(spc.gridPart()),
    spc_(spc),
    mapper_(spc, spc.mapper().containedMapper(), component),
    component_(component),
    baseSetVec_(GeometryIdentifier::numTypes, 0)
  {
    // create info for all geom types 
    AllGeomTypes<typename GridPartType :: IndexSetType, GridType>
      allGeomTypes(spc.gridPart().indexSet());

    // get types for codim 0 
    const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);

    // create base sets for all existing geom types
    for(size_t i=0; i<geomTypes.size(); ++i)
    {
      GeometryIdentifier::IdentifierType id =
                      GeometryIdentifier::fromGeo(geomTypes[i]);

      assert(id >= 0 && id < static_cast<int>(GeometryIdentifier::numTypes));

      if (baseSetVec_[id] == 0) 
      {
        baseSetVec_[id] =
          new BaseFunctionSetType( spc.baseFunctionSet( id ) , component);
      }
    } // end for
  }
  
  template <class CombinedSpaceImp>
  inline
  SubSpace<CombinedSpaceImp>::~SubSpace()
  {
    for(size_t i=0; i<baseSetVec_.size(); ++i)
    {
      delete baseSetVec_[i];
    } // end for
  }
  
  //- class SubBaseFunctionSet
  template <class CombinedSpaceImp>
  template <int diffOrd>
  inline
  void SubBaseFunctionSet<CombinedSpaceImp>::
  evaluate (int baseFunct, 
            const FieldVector<deriType, diffOrd>& diffVariable, 
            const DomainType& x, RangeType& phi ) const 
  {
    // Assumption: dimRange == 1
    bSet_.evaluate(baseFunct, diffVariable, x, tmp_);
    phi[0] = tmp_[component_];
  }
  
  //! evaluate base function at quadrature point
  template <class CombinedSpaceImp>
  template <int diffOrd, class QuadratureType >
  inline
  void SubBaseFunctionSet<CombinedSpaceImp>::
  evaluate (int baseFunct, 
            const FieldVector<deriType, diffOrd> &diffVariable, 
            QuadratureType & quad, 
            int quadPoint, RangeType & phi ) const 
  {
    // Assumption: dimRange == 1
    bSet_.evaluate(baseFunct, diffVariable, quad, quadPoint, tmp_);
    phi[0] = tmp_[component_];
  }
  
  //- class SubMapper
  template <class CombinedSpaceImp>
  inline
  int SubMapper<CombinedSpaceImp>::size() const 
  {
    return mapper_.size();
  }

  template <class CombinedSpaceImp>
  template <class EntityType>
  inline
  int SubMapper<CombinedSpaceImp>::
  mapToGlobal(EntityType& en, int localNum) const
  {
    const int containedGlobal = mapper_.mapToGlobal(en, localNum);

    utilGlobal_.newSize(mapper_.size()); // ok, since pointbased specialisation does nothing for newSize
    return utilGlobal_.combinedDof(containedGlobal, component_);
  }

} // end namespace Dune
