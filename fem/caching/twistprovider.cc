namespace Dune {

  template <class ct, int dim>
  std::map<size_t, std::vector<TwistMapper*> > 
  TwistProvider<ct, dim>::mappers_;

  template <class ct, int dim>
  const int TwistProvider<ct, dim>::offset_ = 5;

  template <class ct, int dim>
  const TwistMapper& 
  TwistProvider<ct, dim>::getTwistMapper(const QuadratureType& quad, int twist)
  {
    IteratorType it = mappers_.find(quad.id());
    if (it == mappers_.end()) {
      it = TwistProvider<ct, dim>::addMapper(quad);
    }
    
    assert((it->second)[twist + offset_]);
    return *(it->second)[twist + offset_];
  }
  
  template <class ct, int dim>
  typename TwistProvider<ct, dim>::IteratorType
  TwistProvider<ct, dim>::addMapper(const QuadratureType& quad) 
  {
    std::auto_ptr<TwistMapperCreator<ct, dim> > creator = 
      TwistProvider<ct, dim>::newCreator(quad);
    
    // The vector of TwistMapper for every possible twist
    std::vector<TwistMapper*> mapperVec(offset_ + creator->maxTwist());
    
    for (int twist = creator->minTwist(); twist < creator->maxTwist(); 
         ++twist) {
      mapperVec[offset_ + twist] = creator->createMapper(twist);
    }

    return mappers_.insert(std::make_pair(quad.id(), mapperVec)).first;
  }
  
  template <class ct, int dim>
  std::auto_ptr<TwistMapperCreator<ct, 1> > 
  TwistProvider<ct, dim>::newCreator(const Quadrature<ct, 1>& quad)
  {
    typedef std::auto_ptr<TwistMapperCreator<ct, 1> > AutoPtrType;
    assert(dim == 1);
   
    return AutoPtrType(new LineTwistMapperCreator<ct>(quad));
  }

  template <class ct, int dim>
  std::auto_ptr<TwistMapperCreator<ct, 2> > 
  TwistProvider<ct, dim>::newCreator(const Quadrature<ct, 2>& quad)
  {
    typedef std::auto_ptr<TwistMapperCreator<ct, 2> > AutoPtrType;
    assert(dim == 2);
    
    switch (quad.geo()) {
      case triangle:
      case simplex:
        return AutoPtrType(new TriangleTwistMapperCreator<ct>(quad));
      case quadrilateral:
      case cube:
        return AutoPtrType(new QuadrilateralTwistMapperCreator<ct>(quad));
      default:
        DUNE_THROW(NotImplemented, 
                   "No creator for given GeometryType exists");   
    } // end switch
    
    assert(false); // should newer get here
    return AutoPtrType(new QuadrilateralTwistMapperCreator<ct>(quad)); // dummy


  }

  template <class ct, int dim>
  const ct TwistMapperCreator<ct, dim>::eps_ = 1.0e-5;

  template <class ct, int dim>
  TwistMapperCreator<ct, dim>::TwistMapperCreator(const QuadratureType& quad,
                                                  int minTwist,
                                                  int maxTwist) :
    quad_(quad),
    mat_(0.),
    minTwist_(minTwist),
    maxTwist_(maxTwist)
  {}

  template <class ct, int dim>
  TwistMapper* TwistMapperCreator<ct, dim>::createMapper(int twist) const 
  {
    TwistMapper* mapper = new TwistMapper();
    mapper->indices_.resize(quad_.nop());
    
    buildTransformationMatrix(twist, mat_);

    for (int i = 0; i < quad_.nop(); ++i) {
      PointType pFace = quad_.point(i);
      CoordinateType c(0.0);
      c[0] = 1.0;
      for (int d = 0; d < dim; ++d) {
        c[0] -= pFace[d];
        c[d+1] = pFace[d];
      }

      PointType pRef(0.);     
      mat_.umtv(c, pRef);

      bool found = false;
      for (int j = 0; j < quad_.nop(); ++j) {
        if (samePoint(pRef, quad_.point(j))) {
          mapper->indices_[i] = j;
          found = true;
          break;
        }
      }
      // * Needs to made more generic for non-symmetric quadratures
      assert(found);
    }
    
    return mapper;
  }

  template <class ct, int dim>
  bool TwistMapperCreator<ct, dim>::samePoint(const PointType& first,
                                              const PointType& second) const
  {
    bool result = true;
    for (int i = 0; i < dim; ++i) {
      result &= (first[i] < (second[i] + eps_) 
                 && first[i] > (second[i] - eps_));
    }
    
    return result;
  }

  template <class ct>
  LineTwistMapperCreator<ct>::
  LineTwistMapperCreator(const QuadratureType& quad) :
    TwistMapperCreator<ct, 1>(quad, 0, 1),
    refElem_(ReferenceElements<ct, dim>::cube(quad.geo()))
  {}

  template <class ct>
  void LineTwistMapperCreator<ct>::
  buildTransformationMatrix(int twist, MatrixType& mat) const 
  {
    assert(twist == 0 || twist == 1);
    
    if (twist == 0) {
      mat[0] = refElem_.position(0, 1);
      mat[1] = refElem_.position(1, 1);
    }
    else {
      mat[0] = refElem_.position(1, 1);
      mat[1] = refElem_.position(0, 1);
    }
  }

  template <class ct>
  TriangleTwistMapperCreator<ct>::
  TriangleTwistMapperCreator(const QuadratureType& quad) :
    TwistMapperCreator<ct, 2>(quad, -3, 3),
    refElem_(ReferenceElements<ct, dim>::simplex(quad.geo()))
  {}

  template <class ct>
  void TriangleTwistMapperCreator<ct>::
  buildTransformationMatrix(int twist, MatrixType& mat) const 
  {
   mat = 0.0;

    for (int idx = 0; idx < dim+1; ++idx) {
      int aluIndex = FaceTopo::dune2aluVertex(idx);
      int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }
  }

  template <class ct>
  QuadrilateralTwistMapperCreator<ct>::
  QuadrilateralTwistMapperCreator(const QuadratureType& quad) :
    TwistMapperCreator<ct, 2>(quad, -4, 4),
    refElem_(ReferenceElements<ct, dim>::cube(quad.geo()))
  {}

  template <class ct>
  void QuadrilateralTwistMapperCreator<ct>::
  buildTransformationMatrix(int twist, MatrixType& mat) const 
  {
    mat = 0.0;

    for (int idx = 0; idx < dim+1; ++idx) {
      int aluIndex = FaceTopo::dune2aluVertex(idx);
      int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }
  }

} // end namespace Dune
