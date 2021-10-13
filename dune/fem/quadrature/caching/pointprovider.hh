#ifndef DUNE_FEM_POINTPROVIDER_HH
#define DUNE_FEM_POINTPROVIDER_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/math.hh>

#include <dune/fem/storage/singleton.hh>

//- Local includes
#include "pointmapper.hh"

namespace Dune
{

  namespace Fem
  {

    template< class ct, int dim, int codim >
    class PointProvider
    {
      static_assert( (codim >= 0) && (codim <= 1),
                          "PointProvider exists only for codimension 0 and 1." );
    };

    template <class ct, int dim>
    class PointProvider<ct, dim, 0>
    {
      typedef PointProvider<ct, dim, 0> ThisType;

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

      // points container holding quadrature points
      PointContainerType points_;

      static PointContainerType& points()
      {
        return instance().points_;
      }

      static ThisType& instance()
      {
        return Singleton< ThisType > :: instance();
      }
    };

    // * Add elemGeo later
    template <class ct, int dim>
    class PointProvider<ct, dim, 1>
    {
      typedef PointProvider<ct, dim, 1> ThisType;

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
      typedef std::pair< MapperVectorType, MapperVectorType >  MapperVectorPairType;

    public:
      inline
      static const MapperVectorPairType& getMappers(const QuadratureType& quad,
                                                    const GeometryType& elementGeo);
      // Access for non-symmetric quadratures
      inline
      static const MapperVectorPairType& getMappers(const QuadratureType& quad,
                                                    const LocalPointVectorType& pts,
                                                    const GeometryType& elementGeo);
      inline
      static const GlobalPointVectorType& getPoints(const size_t id,
                                                    const GeometryType& elementGeo);

    private:
      typedef std::map<const QuadratureKeyType, GlobalPointVectorType> PointContainerType;
      typedef std::map<const QuadratureKeyType, MapperVectorPairType > MapperContainerType;

      typedef typename PointContainerType::iterator PointIteratorType;
      typedef typename MapperContainerType::iterator MapperIteratorType;

    private:
      inline
      static MapperIteratorType addEntry(const QuadratureType& quad,
                                         const LocalPointVectorType& pts,
                                         GeometryType elementGeo);

    private:
      // points container holding quadrature points
      PointContainerType   points_;
      // mapper container holding mapping info
      MapperContainerType  mappers_;

      static PointContainerType& points() { return instance().points_; }
      static MapperContainerType& mappers() { return instance().mappers_; }

      static ThisType& instance()
      {
        return Singleton< ThisType > :: instance();
      }
    };

  } // namespace Fem

} // namespace Dune

#include "pointprovider.cc"

#endif // #ifndef DUNE_FEM_POINTPROVIDER_HH
