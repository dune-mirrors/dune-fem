
namespace Dune {
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  const int CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::spaceId_ = 13;

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  CombinedSpace(ContainedDiscreteFunctionSpaceType& spc) :
    BaseType(spc.gridPart(),spaceId_),
    spc_(spc),
    mapper_(spc_, spc_.mapper()),
    baseSetVec_(GeometryIdentifier::numTypes, 0),
    dm_(DofManagerFactoryType::getDofManager(spc_.grid()))
  {
    // get types for codim 0  
    const std::vector<GeometryType>& geomTypes =
                spc_.indexSet().geomTypes(0) ;

    int maxNumDofs = -1;
    // create mappers and base sets for all existing geom types
    for(size_t i=0; i<geomTypes.size(); ++i)
    {
      GeometryIdentifier::IdentifierType id =
                  GeometryIdentifier::fromGeo(geomTypes[i]);
      if(baseSetVec_[id] == 0 )
      {
        baseSetVec_[id] = new BaseFunctionSetType(spc_.getBaseFunctionSet(id));
        maxNumDofs = std::max(maxNumDofs,baseSetVec_[id]->numBaseFunctions());
      }
    }
  }
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  CombinedSpace(GridPartType& gridpart) :
    BaseType(gridpart), // ,spaceId_),
    spc_(gridpart),
    mapper_(spc_, spc_.mapper()),
    baseSetVec_(GeometryIdentifier::numTypes, 0),
    dm_(DofManagerFactoryType::getDofManager(spc_.grid()))
  {
    // get types for codim 0  
    const std::vector<GeometryType>& geomTypes =
                spc_.indexSet().geomTypes(0) ;

    int maxNumDofs = -1;
    // create mappers and base sets for all existing geom types
    for(size_t i=0; i<geomTypes.size(); ++i)
    {
      GeometryIdentifier::IdentifierType id =
                  GeometryIdentifier::fromGeo(geomTypes[i]);
      if(baseSetVec_[id] == 0 )
      {
        baseSetVec_[id] = new BaseFunctionSetType(spc_.getBaseFunctionSet(id));
        maxNumDofs = std::max(maxNumDofs,baseSetVec_[id]->numBaseFunctions());
      }
    }
  }
  
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  CombinedSpace<DiscreteFunctionSpaceImp, N, policy>::
  ~CombinedSpace() 
  {
    for (unsigned int i = 0; i < baseSetVec_.size(); ++i) {
      delete baseSetVec_[i];
      baseSetVec_[i] = 0;
    }
  }

  //- CombinedBaseFunctionSet
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <int diffOrd>
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluate(int baseFunct, 
           const FieldVector<deriType, diffOrd> &diffVariable,
           const DomainType & x, RangeType & phi) const 
  {
    baseFunctionSet_.evaluate(util_.containedDof(baseFunct), 
                              diffVariable, x, containedResult_);
    expand(baseFunct, containedResult_, phi);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <int diffOrd, class QuadratureType>
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluate(int baseFunct, 
           const FieldVector<deriType, diffOrd> &diffVariable, 
           QuadratureType & quad, 
           int quadPoint, RangeType & phi) const 
  {
    baseFunctionSet_.evaluate(util_.containedDof(baseFunct), diffVariable,
                              quad, quadPoint, containedResult_);
    expand(baseFunct, containedResult_, phi);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  expand(int baseFunct, const ContainedRangeType& arg, RangeType& dest) const 
  {
    dest = 0.0;
    assert(arg.dim() == 1); // only DimRange == 1 allowed
    dest[util_.component(baseFunct)] = arg[0];
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateScalar(int baseFunct, 
                 const DomainType& x, 
                 ContainedRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < baseFunctionSet_.numBaseFunctions());
    baseFunctionSet_.eval(baseFunct, x, phi);
  }
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class QuadratureType> 
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateScalar(int baseFunct, 
                 const QuadratureType & quad, int p, 
                 ContainedRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < baseFunctionSet_.numBaseFunctions());
    baseFunctionSet_.eval(baseFunct, quad, p, phi);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class QuadratureType> 
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  jacobianScalar(int baseFunct, 
                 const QuadratureType & quad, int p, 
                 ContainedJacobianRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < baseFunctionSet_.numBaseFunctions());
    baseFunctionSet_.jacobian(baseFunct, quad,p, phi);
  }
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline void CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  jacobianScalar(int baseFunct, 
                 const DomainType& x, 
                 ContainedJacobianRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < baseFunctionSet_.numBaseFunctions());
    baseFunctionSet_.jacobian(baseFunct, x, phi);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  inline typename CombinedBaseFunctionSet<DiscreteFunctionSpaceImp,N,policy>::DofType
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateSingle(int baseFunct, 
                 const DomainType& xLocal,
                 const RangeType& factor) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numBaseFunctions() );

    baseFunctionSet_.eval(util_.containedDof(baseFunct), xLocal, phi_ );
    assert( util_.component(baseFunct) < RangeType :: dimension );
    return factor[util_.component(baseFunct)]*phi_[0];
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class QuadratureType> 
  inline typename CombinedBaseFunctionSet<DiscreteFunctionSpaceImp,N,policy>::DofType
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateSingle(int baseFunct, 
                 const QuadratureType & quad, int qp, 
                 const RangeType& factor) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numBaseFunctions() );
    baseFunctionSet_.eval(util_.containedDof(baseFunct), quad, qp, phi_);
    assert( util_.component(baseFunct) < RangeType :: dimension );
    return factor[util_.component(baseFunct)]*phi_[0];
    //return evaluateSingle(baseFunct,quad.point(qp),factor);
  }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class Entity>
  inline typename
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::DofType
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateGradientSingle(int baseFunct,
                         Entity& en,
                         const DomainType& xLocal,
                         const JacobianRangeType& factor) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numBaseFunctions() );

    //baseFunctionSet_.jacobian(baseFunct, x, phi);
    gradScaled_ = 0.0;
    
    baseFunctionSet_.jacobian(util_.containedDof(baseFunct), xLocal, grad_ );
    en.geometry().jacobianInverseTransposed(xLocal).
      umv(grad_[0], gradScaled_ );
    //! is this right?
    //return factor[util_.component(baseFunct)]*jTmp[0];
    assert( util_.component(baseFunct) >= 0 );
    assert( util_.component(baseFunct) < JacobianRangeType :: rows );
    
    return gradScaled_ * factor[util_.component(baseFunct)];
 }

  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class Entity, class QuadratureType>
  inline typename
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::DofType
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateGradientSingle(int baseFunct,
                         Entity& en,
                         const QuadratureType & quad, int qp, 
                         const JacobianRangeType& factor) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numBaseFunctions() );

    //baseFunctionSet_.jacobian(baseFunct, x, phi);
    gradScaled_ = 0.0;
    
    baseFunctionSet_.jacobian(util_.containedDof(baseFunct), quad, qp, grad_ );
    en.geometry().jacobianInverseTransposed( quad.point(qp) ).
      umv(grad_[0], gradScaled_ );
    //! is this right?
    //return factor[util_.component(baseFunct)]*jTmp[0];
    assert( util_.component(baseFunct) >= 0 );
    assert( util_.component(baseFunct) < JacobianRangeType :: rows );
    
    return gradScaled_ * factor[util_.component(baseFunct)];
  }
  template <class DiscreteFunctionSpaceImp, int N, DofStoragePolicy policy>
  template <class Entity, class QuadratureType>
  inline typename
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::DofType
  CombinedBaseFunctionSet<DiscreteFunctionSpaceImp, N, policy>::
  evaluateGradientTransformed(int baseFunct,
                              Entity& en,
                              const QuadratureType & quad, int qp, 
                              const JacobianRangeType& factor) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numBaseFunctions() );

    //baseFunctionSet_.jacobian(baseFunct, x, phi);a
    baseFunctionSet_.jacobian(util_.containedDof(baseFunct), quad, qp, grad_ );
    /*
    gradScaled_ = 0.0;
    
    en.geometry().jacobianInverseTransposed( quad.point(qp) ).
      umv(grad_[0], gradScaled_ );
    //! is this right?
    //return factor[util_.component(baseFunct)]*jTmp[0];
    assert( util_.component(baseFunct) >= 0 );
    assert( util_.component(baseFunct) < JacobianRangeType :: rows );
    return gradScaled_ * factor[util_.component(baseFunct)];
    */
    return grad_ [0]* factor[util_.component(baseFunct)];
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
    assert( false );
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
    assert(false);
    return newIndex(num);
    /*
    DofConversionUtility<policy> 
      tmpUtilGlobal(chooseSize(N, mapper_.oldSize(), Int2Type<policy>()));

    const int component = tmpUtilGlobal.component(num);
    const int contained = tmpUtilGlobal.containedDof(num);

    const int containedNew = mapper_.oldIndex(contained);

    return tmpUtilGlobal.combinedDof(containedNew, component);
    */
  }
 
} // end namespace Dune
