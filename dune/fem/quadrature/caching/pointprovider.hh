#ifndef DUNE_FEM_POINTPROVIDER_HH
#define DUNE_FEM_POINTPROVIDER_HH

//- System includes
#include <vector>
#include <map>
#include <unordered_map>

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
      static void registerQuadrature(const QuadratureType& quad)
      {
        instance().registerQuadratureImpl( quad );
      }

      inline
      static const GlobalPointVectorType& getPoints(const size_t id,
                                                    const GeometryType& elementGeo)
      {
        return instance().getPointsImpl( id, elementGeo );
      }

      PointProvider() : threadPool_( MPIManager::threadPool() ), points_() {}

    private:
      PointProvider( const PointProvider& ) = delete;

      inline void registerQuadratureImpl(const QuadratureType& quad);

      inline const GlobalPointVectorType& getPointsImpl(const size_t id,
                                                        const GeometryType& elementGeo);

      typedef std::unordered_map<QuadratureKeyType, GlobalPointVectorType> PointContainerType;
      typedef typename PointContainerType::iterator PointIteratorType;

      // points container holding quadrature points
      const typename MPIManager::ThreadPoolType& threadPool_;
      PointContainerType points_;

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

    private:
      typedef std::map<const QuadratureKeyType, GlobalPointVectorType> PointContainerType;
      typedef std::map<const QuadratureKeyType, MapperVectorPairType > MapperContainerType;

      typedef typename PointContainerType::iterator PointIteratorType;
      typedef typename MapperContainerType::iterator MapperIteratorType;

    public:
      inline
      static const MapperVectorPairType& getMappers(const QuadratureType& quad,
                                                    const GeometryType& elementGeo)
      {
        return instance().getMappersImpl( quad, elementGeo );
      }

      // Access for non-symmetric quadratures
      inline
      static const MapperVectorPairType& getMappers(const QuadratureType& quad,
                                                    const LocalPointVectorType& pts,
                                                    const GeometryType& elementGeo)
      {
        return instance().getMappersImpl( quad, pts, elementGeo );
      }

      inline
      static const GlobalPointVectorType& getPoints(const size_t id,
                                                    const GeometryType& elementGeo)
      {
        return instance().getPointsImpl( id, elementGeo );
      }

      inline
      const GlobalPointVectorType& getPointsImpl(const size_t id,
                                                 const GeometryType& elementGeo);

      inline const MapperVectorPairType& getMappersImpl(const QuadratureType& quad,
                                                        const GeometryType&
                                                        elementGeo);

      inline
      const MapperVectorPairType& getMappersImpl(const QuadratureType& quad,
                                                 const LocalPointVectorType& pts,
                                                 const GeometryType& elementGeo);

    private:
      inline
      static MapperIteratorType addEntry(const QuadratureType& quad,
                                         const LocalPointVectorType& pts,
                                         GeometryType elementGeo)
      {
        instance().addEntryImpl(quad, pts, elementGeo);
      }

      inline
      MapperIteratorType addEntryImpl(const QuadratureType& quad,
                                      const LocalPointVectorType& pts,
                                      GeometryType elementGeo);

    private:
      const typename MPIManager::ThreadPoolType& threadPool_;
      // points container holding quadrature points
      PointContainerType   points_;
      // mapper container holding mapping info
      MapperContainerType  mappers_;

      PointProvider( const PointProvider& ) = delete;
    public:
      // this should only be called from the Singleton class
      PointProvider() : threadPool_( MPIManager::threadPool() ), points_(), mappers_() {}

      static ThisType& instance()
      {
        return Singleton< ThisType > :: instance();
      }
    };

  } // namespace Fem

} // namespace Dune

#include "pointprovider.cc"

#endif // #ifndef DUNE_FEM_POINTPROVIDER_HH
