namespace Dune
{

  namespace Fem
  {

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
    const typename TwistProvider<ct, dim>::TwistStorageType&
    TwistProvider<ct, dim>::getTwistStorage(const QuadratureType& quad)
    {
      typedef const TwistStorageType* TwistStoragePtr;
      assert( quad.id() < MapperContainer::instance().size() );
      TwistStoragePtr& ptr = MapperContainer::instance()[ quad.id() ];

      if( ptr == 0 )
      {
        TwistMapperCreator<ct, dim> creator(quad);
        ptr = creator.createStorage();
      }

      assert( ptr != 0 );
      return *ptr ;
    }

    template <class ct, int dim>
    const ct TwistMapperCreator<ct, dim>::eps_ = 1.0e-5;

    template <class ct, int dim>
    TwistMapperCreator<ct, dim>::TwistMapperCreator(const QuadratureType& quad) :
      quad_(quad),
      helper_( 0 )
    {
      const GeometryType geoType = quad.geometryType();
      if (dim == 0) {
        helper_ = new PointTwistMapperStrategy<ct,dim>( geoType );
      }
      else if (dim == 1) {
        helper_ = new LineTwistMapperStrategy<ct, dim>( geoType );
      }
      else
      {
        assert (dim == 2);

        if(geoType.isTriangle())
        {
          helper_ = new TriangleTwistMapperStrategy<ct, dim>( geoType );
          return ;
        }
        if( geoType.isQuadrilateral())
        {
          helper_ = new QuadrilateralTwistMapperStrategy<ct,dim>( geoType );
          return ;
        }
        DUNE_THROW(NotImplemented,
                   "No creator for given GeometryType exists");
      }
    }

    template <class ct, int dim>
    TwistMapperCreator<ct, dim>::~TwistMapperCreator()
    {
      delete helper_ ; helper_ = 0 ;
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
      for (int twist = helper_->minTwist();
           twist < helper_->maxTwist(); ++twist)
      {
        MapperType mapper(quad_.nop());

        const MatrixType mat = helper_->buildTransformationMatrix(twist);

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

  } // namespace Fem

} // namespace Dune
