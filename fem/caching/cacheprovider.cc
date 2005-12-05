

namespace Dune {

  template <class GridImp, int codim>
  const CacheProvider<GridImp, codim>::CachePointMapperType&
  CacheProvider<GridImp, codim>::getMapper(const CachingQuadratureType& quad) 
  {
    CacheProvider::getMapper(quad, Int2Type<unstructured>())
  }

  template <class GridImp, int codim>
  const CacheProvider<GridImp, codim>::CachePointMapperType&
  CacheProvider<GridImp, codim>::getMapper(const CachingQuadratureType& quad,
                                           Int2Type<true>) 
  {
    
  }

  template <class GridImp, int codim>
  const CacheProvider<GridImp, codim>::CachePointMapperType&
  CacheProvider<GridImp, codim>::getMapper(const CachingQuadratureType& quad,
                                           Int2Type<false>) 
  {
 
    
  }

  

}
