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
    TwistMapperCreator<ct, dim> creator(quad);
    
    // The vector of TwistMapper for every possible twist
    std::vector<TwistMapper*> mapperVec(offset_ + creator.maxTwist());
    
    for (int twist = creator.minTwist(); twist < creator.maxTwist(); 
         ++twist) {
      mapperVec[offset_ + twist] = creator.createMapper(twist);
    }

    return mappers_.insert(std::make_pair(quad.id(), mapperVec)).first;
  }
  
  template <class ct, int dim>
  const ct TwistMapperCreator<ct, dim>::eps_ = 1.0e-5;

  template <class ct, int dim>
  TwistMapperCreator<ct, dim>::TwistMapperCreator(const QuadratureType& quad) :
    quad_(quad),
    helper_()
  {
    typedef std::auto_ptr<TwistMapperStrategy<ct, dim> > AutoPtrType;
    
    if (dim == 1) {
      helper_ = 
        AutoPtrType(new LineTwistMapperStrategy<ct, dim>(quad.geo()));
    } 
    else {
      assert (dim == 2);

      switch (quad.geo()) {
      case triangle:
      case simplex:
        helper_ = 
          AutoPtrType(new TriangleTwistMapperStrategy<ct, dim>(quad.geo()));
      case quadrilateral:
      case cube:
        helper_ = 
         AutoPtrType(new QuadrilateralTwistMapperStrategy<ct,dim>(quad.geo()));
      default:
        DUNE_THROW(NotImplemented, 
                   "No creator for given GeometryType exists");
      } // end switch
    }
  }

  template <class ct, int dim>
  TwistMapper* TwistMapperCreator<ct, dim>::createMapper(int twist) const 
  {
    TwistMapper* mapper = new TwistMapper(quad_.nop());
    
    const MatrixType& mat = helper_->buildTransformationMatrix(twist);

    for (int i = 0; i < quad_.nop(); ++i) {
      PointType pFace = quad_.point(i);
      CoordinateType c(0.0);
      c[0] = 1.0;
      for (int d = 0; d < dim; ++d) {
        c[0] -= pFace[d];
        c[d+1] = pFace[d];
      }

      PointType pRef(0.);     
      mat.umtv(c, pRef);

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

  template <class ct, int dim>
  LineTwistMapperStrategy<ct, dim>::
  LineTwistMapperStrategy(GeometryType geo) :
    TwistMapperStrategy<ct, dim>(0, 1),
    refElem_(ReferenceElements<ct, dim>::cube(geo)),
    mat_(0.)
  {
    assert(dim == 1);
  }

  template <class ct, int dim>
  const typename TwistMapperStrategy<ct, dim>::MatrixType&
  LineTwistMapperStrategy<ct, dim>::buildTransformationMatrix(int twist) const 
  {
    assert(twist == 0 || twist == 1);
    
    if (twist == 0) {
      mat_[0] = refElem_.position(0, 1);
      mat_[1] = refElem_.position(1, 1);
    }
    else {
      mat_[0] = refElem_.position(1, 1);
      mat_[1] = refElem_.position(0, 1);
    }
  }

  template <class ct, int dim>
  TriangleTwistMapperStrategy<ct, dim>::
  TriangleTwistMapperStrategy(GeometryType geo) :
    TwistMapperStrategy<ct, dim>(-3, 3),
    refElem_(ReferenceElements<ct, dim>::simplices(geo)),
    mat_(0.)
  {
    assert(dim == 2);
  }

  template <class ct, int dim>
  const typename TwistMapperStrategy<ct, dim>::MatrixType&
  TriangleTwistMapperStrategy<ct, dim>::
  buildTransformationMatrix(int twist) const 
  {
   mat_ = 0.0;

    for (int idx = 0; idx < dim+1; ++idx) {
      int aluIndex = FaceTopo::dune2aluVertex(idx);
      int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat_[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }
  }

  template <class ct, int dim>
  QuadrilateralTwistMapperStrategy<ct, dim>::
  QuadrilateralTwistMapperStrategy(GeometryType geo) :
    TwistMapperStrategy<ct, dim>(-4, 4),
    refElem_(ReferenceElements<ct, dim>::cube(geo)),
    mat_(0.)
  {
    assert(dim == 2);
  }

  template <class ct, int dim>
  const typename TwistMapperStrategy<ct, dim>::MatrixType& 
  QuadrilateralTwistMapperStrategy<ct, dim>::
  buildTransformationMatrix(int twist) const 
  {
    mat_ = 0.0;

    for (int idx = 0; idx < dim+1; ++idx) {
      int aluIndex = FaceTopo::dune2aluVertex(idx);
      int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat_[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }
  }

} // end namespace Dune
