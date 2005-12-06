namespace Dune {

  template <class ct, int dim>
  typename TwistProvider<ct, dim>::MapperContainerType
  TwistProvider<ct, dim>::mappers_;

  template <class ct, int dim>
  const int TwistStorage<ct, dim>::offset_ = 5;

  template <class ct, int dim>
  TwistStorage<ct, dim>::TwistStorage(int maxTwist) :
    mappers_(offset_ + maxTwist),
    points_()
  {}

  template <class ct, int dim>
  void TwistStorage<ct, dim>::addMapper(const MapperType& mapper,
                                        int twist)
  {
    mappers_[twist + TwistStorage<ct, dim>::offset_] = mapper;
  }

  template <class ct, int dim>
  void TwistStorage<ct, dim>::addPoint(const PointType& point)
  {
    points_.push_back(point);
  }

  template <class ct, int dim>
  const typename TwistStorage<ct, dim>::MapperType& 
  TwistStorage<ct, dim>::getMapper(int twist) const 
  {
    return mappers_[twist + TwistStorage<ct, dim>::offset_];
  }

  template <class ct, int dim>
  const typename TwistStorage<ct, dim>::PointVectorType&
  TwistStorage<ct, dim>::getPoints() const 
  {
    return points_;
  }

  template <class ct, int dim>
  const typename TwistProvider<ct, dim>::TwistStorageType& 
  TwistProvider<ct, dim>::getTwistStorage(const QuadratureType& quad)
  {
    IteratorType it = mappers_.find(quad.id());
    if (it == mappers_.end()) {
      it = TwistProvider<ct, dim>::createMapper(quad);
    }
    
    assert(it->second);
    return *(it->second);
  }
  
  template <class ct, int dim>
  typename TwistProvider<ct, dim>::IteratorType
  TwistProvider<ct, dim>::createMapper(const QuadratureType& quad) 
  {
    TwistMapperCreator<ct, dim> creator(quad);
    const TwistStorageType* storage = creator.createStorage();

    return mappers_.insert(std::make_pair(quad.id(), storage)).first;
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
  const typename TwistMapperCreator<ct, dim>::TwistStorageType* 
  TwistMapperCreator<ct, dim>::createStorage() const 
  {
    TwistStorageType* storage = new TwistStorageType(helper_->maxTwist());

    // Add quadrature points
    for (int i = 0; i < quad_.nop(); ++i) {
      storage->addPoint(quad_.point(i));
    }

    // Loop over all twists
    for (int twist = helper_->minTwist();twist < helper_->maxTwist();++twist) {
      MapperType mapper(quad_.nop());
      
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
        // find equivalent quadrature point
        for (int j = 0; j < quad_.nop(); ++j) {
          if (samePoint(pRef, quad_.point(j))) {
            mapper[i] = j;
            found = true;
            break;
          }
        }
        // add point if it is not one of the quadrature points
        if (!found) {
          storage->addPoint(pRef);
          mapper.push_back(mapper.size());
        }
      } // for all quadPoints
      storage->addMapper(mapper, twist);
    } // for all twists
    
    return storage;
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
