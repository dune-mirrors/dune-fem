

namespace Dune {

  template <>
  const CacheProvider<>::CachePointMapperType&
  CacheProvider<>::getMapper(const CachingQuadratureType& quad) 
  {
    getMapper(quad, Int2Type<unstructured>())
  }

  template <>
  const CacheProvider<>::CachePointMapperType&
  CacheProvider<>::getMapper(const CachingQuadratureType& quad,
                             Int2Type<true>) 
  {
    // storage problem: im kombinierten Fall, was speichere ich hier?
    // (speichere ich ueberhaupt? doppel-vektor hack funktioniert nicht!)
    
  }

  

}
