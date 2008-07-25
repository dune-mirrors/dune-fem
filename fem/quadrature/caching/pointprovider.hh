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

  template <class ct, int dim, int codim>
  class PointProvider 
  {
    typedef CompileTimeChecker<false> 
    Point_Provider_exists_only_for_codims_1_and_2;
  };

  template <class ct, int dim>
  class PointProvider<ct, dim, 0>
  {
  private:
    typedef CachingTraits<ct, dim> Traits;

  public:
    typedef typename Traits::QuadratureType QuadratureType;
    typedef typename Traits::PointVectorType GlobalPointVectorType;
    typedef typename Traits::QuadratureKeyType QuadratureKeyType;
    
  public:
    inline
    static void registerQuadrature(const QuadratureType& quad);
    
    inline
    static const GlobalPointVectorType& getPoints(const size_t id,
                                                  const GeometryType& elementGeo);

  private:
    typedef std::map<const QuadratureKeyType, GlobalPointVectorType> PointContainerType;
    typedef typename PointContainerType::iterator PointIteratorType;

  private:
    static PointContainerType points_;
  };

  // * Add elemGeo later
  template <class ct, int dim>
  class PointProvider<ct, dim, 1>
  {
    enum { codim = 1 };
    typedef CachingTraits<ct, dim-codim> Traits;
 
  public:
    typedef typename Traits::QuadratureType QuadratureType;
    typedef typename Traits::PointType LocalPointType;
    typedef typename Traits::PointVectorType LocalPointVectorType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::MapperVectorType MapperVectorType;
    typedef FieldVector<ct, dim> GlobalPointType;
    typedef std::vector<GlobalPointType> GlobalPointVectorType;
    typedef typename Traits::QuadratureKeyType QuadratureKeyType;
    
  public:
    inline
    static const MapperVectorType& getMappers(const QuadratureType& quad,
                                              const GeometryType& elementGeo);
    // Access for non-symmetric quadratures
    inline
    static const MapperVectorType& getMappers(const QuadratureType& quad,
                                              const LocalPointVectorType& pts,
                                              const GeometryType& elementGeo);
    inline
    static const GlobalPointVectorType& getPoints(const size_t id,
                                                  const GeometryType& elementGeo);
    
  private:
    typedef std::map<const QuadratureKeyType, GlobalPointVectorType> PointContainerType;
    typedef std::map<const QuadratureKeyType, MapperVectorType> MapperContainerType;
    typedef typename PointContainerType::iterator PointIteratorType;
    typedef typename MapperContainerType::iterator MapperIteratorType;

  private:
    inline
    static MapperIteratorType addEntry(const QuadratureType& quad,
                                       const LocalPointVectorType& pts,
                                       GeometryType elementGeo);
    inline
    static bool sameGeometry(GeometryType geo1, GeometryType geo2);

  private:
    static PointContainerType points_;
    static MapperContainerType mappers_;
  };

} // end namespace Dune

#include "pointprovider.cc"

#endif
