#ifndef DUNE_CACHEPROVIDER_HH
#define DUNE_CACHEPROVIDER_HH

#include <vector>

#include "twistprovider.hh"

namespace Dune {

  // * Idee: nur eine pointMapper-Klasse
  // * benuetze dann swap zur initialisierung mit nur einem vector! (-> der wird dann ungueltig)

  template <class ct, int dim>
  class CachePointProvider;

  // * ueberarbeiten...
  template <class GridImp, int codim>
  class CachePointStorage 
  {
  public:
    class Key {
    public:
      Key(int faceIdx, int twist) :
        faceIdx_(faceIdx),
        twist_(twist)
      {}

      bool operator<(const Key& other) const {
        if (faceIdx_ < other.faceIdx_) {
          return true;
        }
        return (twist_ < other.twist_);
      }
      
      int faceIndex() const {
        return faceIndex_;
      }

      int twist() const {
        return twist_;
      }

    private:
      int faceIdx_;
      int twist_;
      // * Need to add that later to deal with multiple element types
      // GeometryType elemGeo_;
    };


  public:

  public:
    CachePointStorage();
    
    void addPoints(const Key& key, const CachePointMapper& points) {
      points_.insert(std::make_pair(key, points));
    }

    const CachePointMapper& getMapper(const Key& key) const {
      assert(points_.find(key) != points.end());
      return points_.find(key);
    }

  private:
    typedef std::map<Key, CachePointMapper> MapType;
    typedef typename MapType::iterator IteratorType;

  private:
    MapType points_;
  };


  template <class GridImp, int codim>
  class CacheProvider {
  public:
    typedef CacheQuadrature<GridImp, codim> CacheQuadratureType;
    typedef CachePointMapper MappingType;

  public:
    static const MappingType& getMapper(const CachingQuadratureType& quad);

  private:
    typedef CachePointStorage<GridImp, codim> CachePointStorageType;

  private:
    static const MappingType& getMapper(const CachingQuadratureType& quad,
                                        Int2Type<true>);
    static const MappingType& getMapper(const CachingQuadratureType& quad,
                                        Int2Type<false>);

  private:
    static std::vector<CachePointStorageType> storage_;
  };

  template <class ct, int dim>
  class CachePointProvider {
  public:
    CachePointProvider(const TwistMapper& twist,
                       const CachePointMapper& points) :
      mapping_(twist.numPoints())
    {
      for (size_t i = 0; i < twist.numPoints(); ++i) {
        mapping_[i] = points.index(twist.index(i));
      }
    }

    CachePointProvider(const CachePointMapper& points) :
      mapping_(points.numPoints())
    {
      for (size_t i = 0; i < points.numPoints(); ++i) {
        mapping_[i] = points.index(i);
      }
    }

    size_t index(size_t quadPoint) const {
      return points_[quadPoint];
    }

  private:
    std::vector<size_t> mapping_;
  };



  // * OLD CODE


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
