#ifndef DUNE_CACHEPROVIDER_HH
#define DUNE_CACHEPROVIDER_HH

#include <vector>

#include "twistprovider.hh"

namespace Dune {

  template <class ct, int dim, bool unstructured>
  class CacheProvider {};
  template <class ct, int dim, bool unstructured>
  class CachePointMapper {};

  // No need for twist correction
  template <class ct, int dim>
  class CacheProvider<ct, dim, false> {
  public:
    typedef Quadrature<ct, dim> QuadratureType;
    typedef CachePointMapper<ct, dim, false> CachePointMapperType;

  public:
    static const CachePointMapperType& getMapper(const QuadratureType& quad,
                                                 int twist,
                                                 GeometryType faceGeo,
                                                 GeometryType elemGeo);
    static const CachePointsType& getPoints();

  private:
  };

  template <class GridImp, int codim>
  class CacheProvider<ct, dim, false> {
  public:
    typedef CacheQuadrature<GridImp, codim> CachingQuadratureType;
    typedef typename CachingQuadratureType::QuadratureType QuadratureType;

  public:
    // * if coordinate trafo is used, this cannot be done anymore
    static const CachePointMapperType& getMapper(const CachingQuadratureType& quad);

    // points sind immer angelegt, da zugehoerige Quadratur angelegt
    static const CachePointsType& getPoints(const CachingQuadratureType& quad);
     
  }

  template <class ct, int dim>
  class CachePointMapper<ct, dim, false> {
  public:
    CachePointMapper();

    size_t index(size_t quadPoint) const {
      return indices_[quadPoint];
    }

  private:
    CachePointMapper(const CachePointMapper&);
    CachePointMapper& operator=(const CachePointMapper&);

  private:
    std::vector<size_t> indices_;
  };

  // Needs twist correction
  template <class ct, int dim>
  class CacheProvider<ct, dim, true> {
  public:
    typedef CachePointMapper<ct, dim, true> CachePoint
  };

  template <class ct, int dim>
  class CachePointMapper<ct, dim, true> {
  public:
    CachePointMapper();

    size_t index(size_t quadPoint) const {
      return indices_[twist_.index(quadPoint)];
    }
x
  private:
    CachePointMapper(const CachePointMapper&);
    CachePointMapper& operator=(const CachePointMapper&);

  private:
    const TwistMapper& twist_;
    std::vector<size_t> indices_;
  }

} // end namespace Dune

#endif
