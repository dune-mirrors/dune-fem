namespace Dune {

  template <class ct, int dim>
  PointProvider<ct, dim, 1>::MapperStorageType
  PointProvider<ct, dim, 1>::mappers_;

  template <class ct, int dim>
  PointProvider<ct, dim, 1>::PointStorageType
  PointProvider<ct, dim, 1>::points_;

  template <class ct, int dim>
  const PointProvider<ct, dim, 1>::PointMapper&
  PointProvider<ct, dim, 1>::getMapper(const QuadratureType& quad, int faceIdx)
  {
    MapperIteratorType it = PointProvider::mappers_.find(quad.id());
    
    if (it == mappers_.end()) {
      it = addEntry();
    }

    return *(it->second);
  }

  template <class ct, int dim>
  const PointProvider<ct, dim, 1>::PointType&
  PointProvider<ct, dim, 1>::getPoints(size_t id)
  {
    assert(PointProvider::points_.find(id) != PointProvider::points_.end());
    // * more to come
  }

  template <class ct, int dim>
  void PointProvider<ct, dim, 1>::
  addEntry(size_t id, GeometryType elemGeo, const Points& points) 
  {
    // * more to come
  }
  
}
