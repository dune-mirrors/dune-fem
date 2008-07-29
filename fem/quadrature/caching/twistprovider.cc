namespace Dune {

  template <class ct, int dim>
  TwistStorage<ct, dim>::TwistStorage(int minTwist, int maxTwist) :
    mappers_(Traits::twistOffset_ + maxTwist),
    points_(),
    minTwist_(minTwist),
    maxTwist_(maxTwist)
  {}

  template <class ct, int dim>
  void TwistStorage<ct, dim>::addMapper(const MapperType& mapper,
                                        int twist)
  {
    mappers_[twist + Traits::twistOffset_] = mapper;
  }

  template <class ct, int dim>
  size_t TwistStorage<ct, dim>::addPoint(const PointType& point)
  {
    size_t indexOfPoint = points_.size();
    points_.push_back(point);

    return indexOfPoint;
  }

  template <class ct, int dim>
  const typename TwistStorage<ct, dim>::MapperType& 
  TwistStorage<ct, dim>::getMapper(int twist) const 
  {
    return mappers_[twist + Traits::twistOffset_];
  }

  template <class ct, int dim>
  const typename TwistStorage<ct, dim>::PointVectorType&
  TwistStorage<ct, dim>::getPoints() const 
  {
    return points_;
  }

  template <class ct, int dim>
  int TwistStorage<ct, dim>::minTwist() const 
  {
    return minTwist_;
  }

  template <class ct, int dim>
  int TwistStorage<ct, dim>::maxTwist() const 
  {
    return maxTwist_;
  }

  template <class ct, int dim>
  const typename TwistProvider<ct, dim>::TwistStorageType& 
  TwistProvider<ct, dim>::getTwistStorage(const QuadratureType& quad)
  {
    MapperContainerType& mappers = MapperContainer::instance();
    IteratorType it = mappers.find(quad.id());
    if (it == mappers.end()) {
      it = TwistProvider<ct, dim>::createMapper(quad);
    }
    
    assert(it->second);
    return *(it->second);
  }
  
  template <class ct, int dim>
  typename TwistProvider<ct, dim>::IteratorType
  TwistProvider<ct, dim>::createMapper(const QuadratureType& quad) 
  {
    MapperContainerType& mappers = MapperContainer::instance();
    TwistMapperCreator<ct, dim> creator(quad);
    const TwistStorageType* storage = creator.createStorage();

    return mappers.insert(std::make_pair(quad.id(), storage)).first;
  }
  
  template <class ct, int dim>
  const ct TwistMapperCreator<ct, dim>::eps_ = 1.0e-5;

  template <class ct, int dim>
  TwistMapperCreator<ct, dim>::TwistMapperCreator(const QuadratureType& quad) :
    quad_(quad),
    helper_()
  {
    typedef std::auto_ptr<TwistMapperStrategy<ct, dim> > AutoPtrType;
    
    if (dim == 0) {
      helper_ = AutoPtrType(new
      PointTwistMapperStrategy<ct,dim>(quad.geometry()));
    }
    else if (dim == 1) {
      helper_ = 
        AutoPtrType(new LineTwistMapperStrategy<ct, dim>(quad.geometry()));
    } 
    else 
    {
      assert (dim == 2);

      const GeometryType geoType = quad.geometry();
      if(geoType.isTriangle()) 
      {
        helper_ = 
          AutoPtrType(new TriangleTwistMapperStrategy<ct, dim>( geoType ) );
        return ;
      }
      if( geoType.isQuadrilateral())
      {
        helper_ = 
         AutoPtrType(new QuadrilateralTwistMapperStrategy<ct,dim>( geoType ) );
        return ;
      }
      DUNE_THROW(NotImplemented, 
                 "No creator for given GeometryType exists");
    }
  }

  template <class ct, int dim>
  const typename TwistMapperCreator<ct, dim>::TwistStorageType* 
  TwistMapperCreator<ct, dim>::createStorage() const 
  {
    TwistStorageType* storage = 
      new TwistStorageType(helper_->minTwist(), helper_->maxTwist());

    // Add quadrature points
    for (size_t i = 0; i < quad_.nop(); ++i) {
      storage->addPoint(quad_.point(i));
    }

    // Loop over all twists
    // and find all mapped quad points 
    for (int twist = helper_->minTwist();twist < helper_->maxTwist();++twist) {
      MapperType mapper(quad_.nop());
      
      const MatrixType& mat = helper_->buildTransformationMatrix(twist);

      for (size_t i = 0; i < quad_.nop(); ++i) 
      {
        // get local quad point 
        PointType pFace = quad_.point(i);

        // create barycentric coordinate
        CoordinateType c(0.0);
        c[0] = 1.0;
        for (int d = 0; d < dim; ++d) {
          c[0] -= pFace[d];
          c[d+1] = pFace[d];
        }

        // apply mapping 
        PointType pRef(0.);     
        mat.umtv(c, pRef);

        bool found = false;
        // find equivalent quadrature point
        for (size_t j = 0; j < quad_.nop(); ++j) 
        {
          // check if | pRef - quad_j | < eps 
          if( (pRef - quad_.point(j)).infinity_norm() < eps_)
          {
            mapper[i] = j;
            found = true;
            break;
          }
        }

        // add point if it is not one of the quadrature points
        if (!found) {
          //if point not found, something is wrong,
          //could be the quadratue or the mapping 
          assert( found ); 
          std::cout << "TwistMapperCreator<ct, dim>::createStorage failed! in: "<<__FILE__<<" line: " << __LINE__ <<"\n";
          abort();
        }
        
      } // for all quadPoints
      storage->addMapper(mapper, twist);
    } // for all twists
    
    return storage;
  }

  template <class ct, int dim>
  PointTwistMapperStrategy<ct, dim>::
  PointTwistMapperStrategy(GeometryType geo) :
    TwistMapperStrategy<ct, dim>(0, 1),
    refElem_(ReferenceElements<ct, dim>::cube(geo)),
    mat_(0.)
  {
    assert(dim == 0);
  }

  template <class ct, int dim>
  const typename TwistMapperStrategy<ct, dim>::MatrixType&
  PointTwistMapperStrategy<ct, dim>::buildTransformationMatrix(int twist) const 
  {
    assert(twist == 0);
    mat_[0] = refElem_.position(0, 1);
    return mat_;
  }
  template <class ct, int dim>
  LineTwistMapperStrategy<ct, dim>::
  LineTwistMapperStrategy(GeometryType geo) :
    TwistMapperStrategy<ct, dim>(0, 2),
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
    return mat_;
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

    for (int idx = 0; idx < dim+1; ++idx) 
    {
      const int aluIndex = FaceTopo::dune2aluVertex(idx);
      const int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat_[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }
    
    return mat_;
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
      const int aluIndex = FaceTopo::dune2aluVertex(idx);
      const int twistedDuneIndex = FaceTopo::alu2duneVertex(aluIndex, twist);
      mat_[idx] = refElem_.position(twistedDuneIndex, dim); // dim == codim here
    }

    return mat_;
  }

} // end namespace Dune
