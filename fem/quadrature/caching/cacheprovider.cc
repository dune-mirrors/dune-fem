namespace Dune {

  template <class GridImp>
  typename CacheProvider<GridImp, 1>::MapperContainerType
  CacheProvider<GridImp, 1>::mappers_;

  template <class GridImp>
  typename CacheProvider<GridImp, 1>::MapperIteratorType
  CacheProvider<GridImp, 1>::createMapper(const QuadratureType& quad,
                                          GeometryType elementGeometry,
                                          Int2Type<true>)
  {
    typedef TwistProvider<ct, dim-codim> TwistProviderType;
    typedef typename TwistProviderType::TwistStorageType TwistStorageType;

    const TwistStorageType& twistMappers =
      TwistProviderType::getTwistStorage(quad);
    const MapperVectorType pointMappers =
      PointProvider<ct, dim, codim>::getMappers(quad, 
                                                twistMappers.getPoints(),
                                                elementGeometry);

    const int numFaces = pointMappers.size();
    const int maxTwist = twistMappers.maxTwist();
    const int minTwist = twistMappers.minTwist();

    QuadratureKeyType key ( elementGeometry, quad.id() );
    MapperIteratorType it = mappers_.insert
      (std::make_pair( key,
                       CacheStorageType(numFaces, maxTwist))).first;

    for (int face = 0; face < numFaces; ++face) 
    {
      for (int twist = minTwist; twist < maxTwist; ++twist) {
        it->second.addMapper(pointMappers[face],
                             twistMappers.getMapper(twist),
                             face, twist);
      }
    }

    return it;
  }

  template <class GridImp>
  typename CacheProvider<GridImp, 1>::MapperIteratorType
  CacheProvider<GridImp, 1>::createMapper(const QuadratureType& quad,
                                          GeometryType elementGeometry,
                                          Int2Type<false>)
  {
    const MapperVectorType pointMappers =
      PointProvider<ct, dim, codim>::getMappers(quad, elementGeometry);

    const int numFaces = pointMappers.size();

    QuadratureKeyType key ( elementGeometry, quad.id() );
    
    MapperIteratorType it = 
      mappers_.insert(std::make_pair(key,
                                     CacheStorageType(numFaces))).first;

    for (int face = 0; face < numFaces; ++face) {
      it->second.addMapper(pointMappers[face], 
                           face);
    }

    return it;
  }

} // end namespace Dune
