
namespace Dune {

  //- AdaptiveDiscreteFunction
  // nothing here, everything in adaptivefunction.hh
  
  
  //- AdaptiveLocalFunction
  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                        DofStorageType& dofVec) :
    spc_(spc),
    dofVec_(dofVec),
    values_(),
    numDofs_(0),
    tmp_(0),
    tmpGrad_(0),
//    tmp_(0.0),  
//    tmpGrad_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}
  
  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  AdaptiveLocalFunction(const ThisType& other) :
    spc_(other.spc_),
    dofVec_(other.dofVec_),
    values_(),
    numDofs_(0),
    tmp_(0.0),
    tmpGrad_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  ~AdaptiveLocalFunction() {}  

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::DofType&
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  operator[] (const int num) 
  {
    assert(init_);
    assert(num >= 0 && num < numDofs());
    // check that storage (dofVec_) and mapper are in sync:
    assert(dofVec_.size() >= spc_.size());
    return (* (values_[num]));
  }
  
  template <class DiscreteFunctionSpaceImp>
  const typename AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::DofType&
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  operator[] (const int num) const 
  {
    assert(init_);
    assert(num >= 0 && num < numDofs());
    // check that storage (dofVec_) and mapper are in sync:
    if( dofVec_.size() != spc_.size() ) 
    assert(dofVec_.size() >= spc_.size());
    return (* (values_[num]));
  }

  template <class DiscreteFunctionSpaceImp>
  int AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  numDofs() const
  {
    assert( numDofs_ == values_.size() );
    return numDofs_;
  }

  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  evaluate(const DomainType& x, RangeType& ret) const 
  {
    assert(init_);
    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      bSet.eval(i, x, tmp_);
      tmp_ *= (*values_[i]);
      ret += tmp_;
    }
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  evaluate(const QuadratureType& quad, 
           const int quadPoint, 
           RangeType& ret) const 
  {
    assert(init_);
    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      bSet.eval(i, quad,quadPoint, tmp_);
      tmp_ *= (*values_[i]);
      ret  += tmp_;
    }
  }

  #if OLDFEM
  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobianLocal(EntityType& en, 
		const DomainType& x, 
		JacobianRangeType& ret) const {
    jacobian(x,ret);
  }

  #endif
  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(EntityType& en, 
	   const DomainType& x, 
	   JacobianRangeType& ret) const
  {
    jacobian(x,ret);
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(EntityType& en, 
           QuadratureType& quad, 
           int quadPoint, 
           JacobianRangeType& ret) const
  {
    jacobian(quad,quadPoint,ret);
  }

  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    assert(init_);
    enum { dim = EntityType::dimension };
    typedef typename DiscreteFunctionSpaceImp::GridType::ctype ctype;

    // get jacobian inverse 
    typedef FieldMatrix<ctype, dim, dim> JacobianInverseType;
    const JacobianInverseType& jti = 
          en().geometry().jacobianInverseTransposed(x);

    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    
    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      // evaluate gradient on reference element
      bSet.jacobian(i, x, tmpGrad_);

      // apply element specific values 
      for (int l = 0; l < dimRange; ++l) 
      {
        tmpGrad_[l] *= *values_[i];
        jti.umv(tmpGrad_[l], ret[l]);
      }
    }    
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(const QuadratureType& quad, 
           const int quadPoint, 
           JacobianRangeType& ret) const
  {
    assert(init_);
    enum { dim = EntityType::dimension };

    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    typedef typename DiscreteFunctionSpaceImp::GridType::ctype ctype;
    
    typedef FieldMatrix<ctype, dim, dim> JacobianInverseType;
    const JacobianInverseType& jti = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

    //JacobianRangeType tmp(0.0);

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      // evaluate gradient on reference element
      bSet.jacobian(i, quad,quadPoint, tmpGrad_);

      // apply element specific values 
      for (int l = 0; l < dimRange; ++l) 
      {
        tmpGrad_[l] *= *(values_[i]);
        jti.umv(tmpGrad_[l], ret[l]);
      }
    }

    /*
    for (int i = 0; i < this->numDofs(); ++i) 
    {
      bSet.jacobian(i, quad,quadPoint, tmpGrad_);
      tmpGrad_ *= *values_[i];
      tmp += tmpGrad_;
    }    
    for (int l = 0; l < dimRange; ++l) 
      jti.umv(tmp[l],ret[l]);
    */
  }

  template <class DiscreteFunctionSpaceImp>
  const typename
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::BaseFunctionSetType& 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::baseFunctionSet() const {
    assert(baseSet_);
    return *baseSet_;
  }
  
  template <class DiscreteFunctionSpaceImp>
  const typename
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::EntityType & 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::en() const 
  {
    assert( en_ );
    return *en_;
  }
  

  // --init
  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  init(const EntityType& en) 
  {
    // NOTE: if init is false, then LocalFunction has been create before. 
    // if multipleGeometryTypes_ is true, then grid has elements 
    // of different geometry type (hybrid grid) and we have to check geometry
    // type again, if not we skip this part, because calling the entity's
    // geometry method is not a cheep call 
    
    if( !init_ || multipleGeometryTypes_ )
    {
      if( geoType_ != en.geometry().type() )
      {
        baseSet_ = &spc_.baseFunctionSet(en);

        numDofs_ = baseSet_->numBaseFunctions();
        values_.resize(numDofs_);

        init_ = true;
        geoType_ = en.geometry().type();
      }
    }

    // cache entity
    en_ = &en;

    assert( geoType_ == en.geometry().type() );
    const int numOfDof = numDofs(); 
    for (int i = 0; i < numOfDof; ++i) 
    {
      values_[i] = &(this->dofVec_[spc_.mapToGlobal(en, i)]);
    }

    return ;
  }
  //- AdaptiveDiscreteFunction (specialisation)
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveDiscreteFunction<
     CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveDiscreteFunction() 
  {
    for (unsigned int i = 0; i < subSpaces_.size(); ++i) {
      delete subSpaces_[i];
      subSpaces_[i] = 0;
    }
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  typename AdaptiveDiscreteFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::SubDiscreteFunctionType
  AdaptiveDiscreteFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  subFunction(int component) 
  {
    SubSpaceType* subSpace = new SubSpaceType(this->space(), component);
    subSpaces_.push_back(subSpace);
    

    return SubDiscreteFunctionType(std::string("Subfunction of ")+this->name(),
                                   *subSpace,
                                   this->dofStorage());
  }

  //- AdaptiveLocalFunction (Specialisation for CombinedSpace)
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                        DofStorageType& dofVec) :
    spc_(spc),
    dofVec_(dofVec),
    values_(),
    numDofs_(0),
    cTmp_(0.0),
    cTmpGradRef_(0.0),
    cTmpGradReal_(0.0),
    tmp_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  AdaptiveLocalFunction(const ThisType& other) :
    spc_(other.spc_),
    dofVec_(other.dofVec_),
    values_(),
    numDofs_(0),
    cTmp_(0.0),
    cTmpGradRef_(0.0),
    cTmpGradReal_(0.0),
    tmp_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveLocalFunction() {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  typename AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::DofType&
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  operator[] (const int num) 
  {
    assert(num >= 0 && num < numDofs());
    return *values_[num/N][static_cast<SizeType>(num%N)];
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::DofType&
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  operator[] (const int num) const 
  {
    assert(num >= 0 && num < numDofs());
    return *values_[num/N][static_cast<SizeType>(num%N)];
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  int AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  numDofs() const 
  {
    //return values_.size()*N;
    assert( numDofs_ == (values_.size()*N) );
    return numDofs_;
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  evaluate(const DomainType& x, 
           RangeType& result) const
  {
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    result = 0.0;
    assert((values_.size()) == bSet.numDifferentBaseFunctions());
    const int valSize = values_.size();
    for (int i = 0; i < valSize; ++i) 
    {
      // Assumption: scalar contained base functions
      bSet.evaluateScalar(i, x, cTmp_);
      for (SizeType j = 0; j < N; ++j) 
      {
        result[j] += cTmp_[0]*(*values_[i][j]);
      }
    }
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template <class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  evaluate(const QuadratureType& quad, 
           const int quadPoint, 
           RangeType& ret) const
  {
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    ret = 0.0;

    assert((values_.size()) == bSet.numDifferentBaseFunctions());
    const int valSize = values_.size();
    for (int i = 0; i < valSize; ++i) 
    {
      // Assumption: scalar contained base functions
      // bSet.evaluateScalar(i, x, cTmp_);
      bSet.evaluateScalar(i, quad,quadPoint, cTmp_);
      for (SizeType j = 0; j < N; ++j) {
        ret[j] += cTmp_[0]*(*values_[i][j]);
      }
    }
  }
  #if OLDFEM
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobianLocal(EntityType& en, 
                const DomainType& x, 
                JacobianRangeType& result) const {
    jacobian(en,x,result);
  }
  #endif
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(EntityType& en, 
	   const DomainType& x, 
	   JacobianRangeType& result) const
  {
    jacobian(x,result);
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(const DomainType& x, JacobianRangeType& result) const
  {
    enum { dim = EntityType::dimension };
    typedef FieldMatrix<DofType, dim, dim> JacobianInverseType;
    result = 0.0;

    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    const JacobianInverseType& jInv = 
      en().geometry().jacobianInverseTransposed(x);

    const int numDiffBaseFct = bSet.numDifferentBaseFunctions();
    for (int i = 0; i < numDiffBaseFct; ++i) 
    {
      //cTmpGradRef_ = 0.0;
      cTmpGradReal_ = 0.0;
      bSet.jacobianScalar(i, x, cTmpGradRef_);
      jInv.umv(cTmpGradRef_[0], cTmpGradReal_[0]);

      for (SizeType j = 0; j < N; ++j) {
        // Assumption: ContainedDimRange == 1
        //cTmpGrad_[0] *= *values_[i][j];
        result[j].axpy(*values_[i][j], cTmpGradReal_[0]);
      }
    }

  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template<class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(EntityType& en, 
           QuadratureType& quad, 
           int quadPoint, 
           JacobianRangeType& ret) const 
  {
    jacobian(quad.point(quadPoint), ret);
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template<class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(const QuadratureType& quad, 
           const int quadPoint, 
           JacobianRangeType& ret) const 
  {
    jacobian(quad.point(quadPoint), ret);
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  assign(int dofNum, const RangeType& dofs) {
    for (SizeType i = 0; i < N; ++i) {
      // Assumption: the local ordering is point based
      *values_[dofNum][i] = dofs[i];
    }
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  int AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  numDifferentBaseFunctions() const 
  {
    return values_.size();
  }

  // --init
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  init(const EntityType& en) 
  {

    // NOTE: if init is false, then LocalFunction has been create before. 
    // if fSpace_.multipleGeometryTypes() is true, then grid has elements 
    // of different geometry type (hybrid grid) and we have to check geometry
    // type again, if not we skip this part, because calling the entity's
    // geometry method is not a cheep call 
    
    if( !init_ || multipleGeometryTypes_ )
    {
      if( geoType_ != en.geometry().type() )
      {
        baseSet_ = &spc_.baseFunctionSet(en);

        numDofs_ = baseSet_->numBaseFunctions();
        values_.resize(numDofs_);

        // real dof number is larger 
        numDofs_ *= N;

        init_ = true;
        geoType_ = en.geometry().type();
      }
    }

    assert( geoType_ == en.geometry().type() );

    // cache entity
    en_ = &en;
    
    const int numDof = numDofs();
    assert( values_.size() == numDof );
    for (int i = 0; i < numDof; ++i) 
    {
      for (SizeType j = 0; j < N; ++j) 
      {
        values_[i][j] = &(dofVec_[spc_.mapToGlobal(en, i*N+j)]);
      } // end for j
    } // end for i
  }
  
  // --baseFunctionSet 
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >:: BaseFunctionSetType& 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  baseFunctionSet() const 
  {
    assert( baseSet_ );
    return *baseSet_;
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >:: EntityType& 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  en() const 
  {
    assert( en_ );
    return *en_;
  }
  
} // end namespace Dune
