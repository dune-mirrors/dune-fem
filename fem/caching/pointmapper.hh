#ifndef DUNE_POINTMAPPER_HH
#define DUNE_POINTMAPPER_HH

#include <vector>

namespace Dune {

  template <class ct, int dim>
  struct CachingTraits {
    typedef Quadrature<ct, dim> QuadratureType;
    typedef typename QuadratureType::CoordinateType PointType;
    typedef std::vector<PointType> PointVectorType;
    typedef std::vector<size_t> MapperType;
    typedef std::vector<MapperType> MapperVectorType;

    static const int twistOffset_ = 5;
  };

}

#endif
