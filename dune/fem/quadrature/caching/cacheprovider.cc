namespace Dune
{
  namespace Fem
  {

    template <class GridPart>
    typename CacheProvider<GridPart, 1>::MapperIteratorType
    CacheProvider<GridPart, 1>::createMapper(const QuadratureType& quad,
                                             GeometryType elementGeometry,
                                             std::integral_constant< bool, true > )
    {
      // make sure we work in single thread mode
      if( ! Fem :: MPIManager :: singleThreadMode() )
      {
        DUNE_THROW(SingleThreadModeError, "CacheProvider::createMapper: only call in single thread mode!");
      }

      typedef TwistProvider<ct, dim-codim> TwistProviderType;
      typedef typename TwistProviderType::TwistStorageType TwistStorageType;

      const TwistStorageType& twistMappers =
        TwistProviderType::getTwistStorage(quad);
      const auto pointMappers =
        PointProvider<ct, dim, codim>::getMappers(quad,
                                                  twistMappers.getPoints(),
                                                  elementGeometry);

      const int numFaces = pointMappers.first.size();
      const int maxTwist = twistMappers.maxTwist();
      const int minTwist = twistMappers.minTwist();

      QuadratureKeyType key ( elementGeometry, quad.id() );

      MapperContainerType& mappers_ = mappers();

      MapperIteratorType it = mappers_.insert
        (std::make_pair( key,
                         CacheStorageType(numFaces, maxTwist))).first;

      for (int face = 0; face < numFaces; ++face)
      {
        for (int twist = minTwist; twist < maxTwist; ++twist)
        {
          it->second.addMapper(pointMappers.first[face],
                               pointMappers.second[face],
                               twistMappers.getMapper(twist),
                               face, twist);
        }
      }

      return it;
    }



    template <class GridPart>
    typename CacheProvider<GridPart, 1>::MapperIteratorType
    CacheProvider<GridPart, 1>::createMapper(const QuadratureType& quad,
                                             GeometryType elementGeometry,
                                             std::integral_constant< bool, false > )
    {
      // make sure we work in single thread mode
      if( ! Fem :: MPIManager :: singleThreadMode() )
      {
        DUNE_THROW(SingleThreadModeError, "CacheProvider::createMapper: only call in single thread mode!");
      }

      const auto pointMappers =
        PointProvider<ct, dim, codim>::getMappers(quad, elementGeometry);

      const int numFaces = pointMappers.first.size();

      QuadratureKeyType key ( elementGeometry, quad.id() );

      MapperContainerType& mappers_ = mappers();

      MapperIteratorType it
        = mappers_.insert(std::make_pair(key, CacheStorageType(numFaces))).first;

      for (int face = 0; face < numFaces; ++face)
        it->second.addMapper(pointMappers.first[face],
                             pointMappers.second[face],
                             face);

      return it;
    }

  } // namespace Fem

} // namespace Dune
