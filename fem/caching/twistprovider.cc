namespace Dune {

  template <class ct, int dim>
  std::map<size_t, TwistMapper*> TwistProvider<ct, dim>::mappers_();

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
    
    assert(*(it->second)[twist + offset_]);
    return *(it->second)[twist + offset_];
  }
  
  template <class ct, int dim>
  typename TwistProvider<ct, dim>::IteratorType
  TwistProvider<ct, dim>::addMapper(const QuadratureType& quad) 
  {
    TwistMapperCreator<ct, dim> creator(quad);
    
    // The vector of TwistMapper for every possible twist
    std::vector<TwistMapper*> mapperVec(offset_ + creator.maxTwist(), 0);
    
    for (int twist = creator.minTwist(); twist < creator.maxTwist(); ++twist) {
      mapperVec_[offset_ + twist] = creator.createMapper(twist);
    }

    return mappers_.insert(quad.id(), mapperVec);
  }
  
  template <class ct, int dim>
  const ct TwistMapperCreator<ct, dim>::eps_ = 1.0e-5;

  template <class ct, int dim>
  TwistMapperCreator<ct, dim>::TwistMapperCreator(const QuadratureType& quad) :
    quad_(quad),
    mat_(0.),
    minTwist_(-1),
    maxTwist_(-1)
  {
    // Recall corners of reference element (at least the needed ones)
    const ReferenceElement& refElem =
      ReferenceElements<ct, dim>::general(quad.geo());
 
    for (int i = 0; i < dim+1; ++i) {
      corners_.push_back(refElem.position(i, dim));
    }

    // Determine minimal and maximal possible twist
    switch (quad.geo()) {
    case line:
      minTwist_ = 0;
      maxTwist_ = 2;
      break;
    case triangle:
      minTwist_ = -3; 
      maxTwist_ = 3;
      break;
    case quadrilateral:
      minTwist_ = -4;
      maxTwist_ = 4;
      break;
    case cube:
      if (dim == 1) {
        minTwist_ = 0;
        maxTwist_ = 2; 
      } 
      else {
        minTwist_ = -4;
        maxTwist_ = 4;   
      }
      break;
    case simplex:
      if (dim == 1) {
        minTwist_ = 0;
        maxTwist_ = 2; 
      }
      else {
        minTwist_ = -3;
        maxTwist_ = 3;
      }
      break;
    default:
      DUNE_THROW(NotImplemented, 
                 "Number of possible twists unknown for geo");
    }
  }

  template <class ct, int dim>
  TwistMapper* TwistMapperCreator<ct, dim>::creatorMapper(int twist) const 
  {
    TwistMapper* mapper = new TwistMapper();
    
    buildTransformationMatrix(twist);

    for (int i = 0; i < quad.nop(); ++i) {
      PointType pRef(0.);
      
      PointType pFace = quad.point(i);
      Coordinate c(0.0);
      c[0] = 1.0;
      for (int d = 0; d < dim; ++d) {
        c[0] -= pFace[d];
        c[d+1] = pFace[d];
      }
     
      mat_.umtv(c, pRef);

      bool found = false;
      for (int j = 0; j < quad.nop(); ++j) {
        if (samePoint(pRef, quad.point(j))) {
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
  void TwistMapperCreator<ct, dim>::buildTransformationMatrix(int twist) const
  {
    mat_ = 0.0;

    for (int idx = 0; idx < dim+1; ++idx) {
      // * Need to consider the different ALU, Dune reference elements
      // * here, the distinction between geometryTypes would come in handy
      int twistedIndex = twist(idx, twist);

      mat_[i] = corners_[twistedIndex];
    }
  }

  template <class ct, int dim>
  bool TwistMapperCreator<ct, dim>::samePoint(const Point& first,
                                              const Point& second) const
  {
    bool result = true;
    for (int i = 0; i < dim; ++i) {
      result &= (first[i] < (second[i] + eps_) 
                 && first[i] > (second[i] - eps_));
    }
    
    return result;
  }

} // end namespace Dune
