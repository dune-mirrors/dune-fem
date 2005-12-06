#ifndef DUNE_POINTPROVIDER_HH
#define DUNE_POINTPROVIDER_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "pointmapper.hh"

namespace Dune {

  template <class ct, int dim>
  class MapperStorage 
  {
  public:
    MapperStorage() {}

    MapperStorage(int numFaces) :
      mappers_(numFaces)
    {}

    void addMapper(const PointMapper& mapper, int face) {
      assert(face >= 0 && face < mappers_.size());
      mappers_[face] = mapper;
    }
    
    const PointMapper& getMapper(int face) const {
      assert(face >= 0 && face < mappers_.size());
      return mappers_[face];
    }
    
  private:
    typedef std::vector<PointMapper> MapperVectorType;

  private:
    MapperVectorType mappers_;
  };

  template <class ct, int dim, int codim>
  class PointProvider 
  {
    typedef CompileTimeChecker<false> 
    Point_Provider_exists_only_for_codims_1_and_2;
  };

  template <class ct, int dim>
  class PointProvider<ct, dim, 0>
  {
  public:
    typedef Quadrature<ct, dim> QuadratureType;
    typedef typename QuadratureType::CoordinateType LocalPointType;
    typedef std::vector<LocalPointType> LocalPointVectorType;
    typedef LocalPointType GlobalPointType;
    typedef LocalPointVectorType GlobalPointVectorType;
    
  public:
    static void addPoints(const QuadratureType& quad, GeometryType elementGeo);

    static const GlobalPointVectorType& getPoints(size_t id,
                                                  GeometryType elementGeo);

  private:
    typedef std::map<size_t, GlobalPointVectorType> PointContainerType;
    typedef typename PointContainerType::iterator PointIteratorType;

  private:
    static PointContainerType points_;
  };

  // * Add elemGeo later
  template <class ct, int dim>
  class PointProvider<ct, dim, 1>
  {
    enum { codim = 1 };
    
  public:
    typedef Quadrature<ct, dim-codim> QuadratureType;
    typedef typename QuadratureType::CoordinateType LocalPointType;
    typedef std::vector<LocalPointType> LocalPointVectorType;
    typedef FieldVector<ct, dim> GlobalPointType;
    typedef std::vector<GlobelPointType> GlobalPointVectorType;
    typedef std::vector<PointMapper> MapperVectorType;
    
  public:
    static const MapperVectorType& getMappers(const QuadratureType& quad,
                                              GeometryType elementGeo);
    // Access for non-symmetric quadratures
    static const MapperVectorType& getMappers(const QuadratureType& quad,
                                              const LocalPointVectorType& pts,
                                              GeometryType elementGeo);

    static const PointVectorType& getPoints(size_t id,
                                            GeometryType elementGeo);
    
  private:
    typedef std::map<size_t, GlobalPointVectorType> PointContainerType;
    typedef std::map<size_t, MapperStorageType> MapperContainerType;
    typedef typename PointContainerType::iterator PointIteratorType;
    typedef typename MapperContainerType::iterator MapperIteratorType;

  private:
    static MapperIteratorType addEntry(const QuadratureType& quad,
                                       const LocalPointVectorType& pts);

  private:
    static PointContainerType points_;
    static MapperContainerType mappers_;
  };

} // end namespace Dune

#include "pointprovider.cc"

#endif
