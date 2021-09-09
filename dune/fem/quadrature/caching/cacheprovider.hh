#ifndef DUNE_FEM_CACHEPROVIDER_HH
#define DUNE_FEM_CACHEPROVIDER_HH

#include <vector>
#include <map>
#include <type_traits>

#include <dune/common/math.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/storage/singleton.hh>

#include "pointmapper.hh"
#include "twistprovider.hh"
#include "pointprovider.hh"

namespace Dune
{

  namespace Fem
  {

    //! Storage class for mappers.
    template <class ct, int dim, bool hasTwists>
    class CacheStorage;



    //! Specialisation for grids with twist (i.e. unstructured ones).
    //---------------------------------------------------------------

    template <class ct, int dim>
    class CacheStorage<ct, dim, true>
    {
    private:
      typedef CachingTraits<ct, dim> Traits;

    public:
      typedef typename Traits::MapperType MapperType;
      typedef typename Traits::MapperPairType MapperPairType;

    public:
      CacheStorage(int numFaces, int maxTwist) :
        mappers_(numFaces)
      {
        for (MapperIteratorType it = mappers_.begin();
             it != mappers_.end(); ++it)
        {
          it->resize(maxTwist + Traits::twistOffset_);
        }
      }

      CacheStorage(const CacheStorage& other) :
        mappers_(other.mappers_)
      {}

      void addMapper(const MapperType& faceMapper,
                     const MapperType& interpolMapper,
                     const MapperType& twistMapper,
                     int faceIndex, int faceTwist)
      {
        assert(twistMapper.size() == faceMapper.size());

        MapperPairType& mapper =
          mappers_[faceIndex][faceTwist + Traits::twistOffset_];
        const size_t size = twistMapper.size();
        mapper.first.resize( size );
        for (size_t i = 0; i < size; ++i)
        {
          mapper.first[i] = faceMapper[twistMapper[i]];
        }

        if( !interpolMapper.empty() )
        {
          assert(twistMapper.size() == interpolMapper.size());
          mapper.second.resize( twistMapper.size() );

          for (size_t i = 0; i < size; ++i) {
            mapper.second[i] = interpolMapper[twistMapper[i]];
          }
        }
      }

      const MapperPairType& getMapper(int faceIndex, int faceTwist) const
      {
        assert( faceTwist + Traits::twistOffset_ >= 0 );
        return mappers_[faceIndex][faceTwist + Traits::twistOffset_];
      }

    private:
      typedef std::vector< std::vector< MapperPairType > > MapperContainerType;
      typedef typename MapperContainerType::iterator MapperIteratorType;

    private:
      MapperContainerType mappers_;
    };



    //! Specialisation for grids without any twists (i.e. Cartesian ones).
    //-------------------------------------------------------------------

    template <class ct, int dim>
    class CacheStorage<ct, dim, false>
    {
    private:
      typedef CachingTraits<ct, dim> Traits;

    public:
      typedef typename Traits::MapperType MapperType;
      typedef typename Traits::MapperPairType MapperPairType;

    public:
      explicit CacheStorage ( int numFaces )
      : mappers_( numFaces )
      {}

      CacheStorage ( const CacheStorage &other )
      : mappers_( other.mappers_ )
      {}

      void addMapper ( const MapperType &mapper,
                       const MapperType &interpolMapper,
                       int faceIndex )
      {
        assert( (faceIndex >= 0) && (faceIndex < (int)mappers_.size()) );
        mappers_[ faceIndex ].first = mapper;
        if( !interpolMapper.empty() )
        {
          assert( interpolMapper.size() == mapper.size() );
          mappers_[ faceIndex ].second = interpolMapper;
        }
      }

      const MapperPairType &getMapper ( int faceIndex, int faceTwist ) const
      {
#ifndef NDEBUG
        if( faceIndex >= (int)mappers_.size() )
          std::cerr << "Error: faceIndex = " << faceIndex << " >= " << mappers_.size() << " = mappers_.size()" << std::endl;
#endif
        assert( (faceIndex >= 0) && (faceIndex < (int)mappers_.size()) );
        return mappers_[ faceIndex ];
      }

    private:
      typedef typename std::vector< MapperPairType > MapperContainerType;

    private:
      MapperContainerType mappers_;
    };


    // CacheProvider
    // -------------

    template< class GridPart, int codim >
    class CacheProvider;

    template <class GridPart>
    class CacheProvider<GridPart, 0>
    {
    private:
      static const int codim = 0;
      static const int dim = GridPart::dimension;
      typedef typename GridPart::ctype ct;
      typedef CachingTraits<ct, dim> Traits;

    public:
      typedef typename Traits::QuadratureType QuadratureType;

    public:
      template <class QuadratureImpl>
      static void registerQuadrature(const QuadratureImpl& quad) {
        // get quadrature implementation
        PointProvider<ct, dim, codim>::registerQuadrature(quad.ipList());
      }
    };

    template <class GridPart>
    class CacheProvider<GridPart, 1>
    {
      typedef CacheProvider<GridPart, 1> ThisType;

      static const int codim = 1;
      static const int dim = GridPart::dimension;
      typedef typename GridPart::ctype ct;
      typedef CachingTraits<ct, dim-codim> Traits;

      // true if grid could have twists
      static const bool hasTwists = ! Dune::Fem::GridPartCapabilities::isCartesian<GridPart>::v ;
    public:
      typedef typename Traits::QuadratureType QuadratureType;
      typedef typename Traits::MapperType MapperType;
      typedef typename Traits::QuadratureKeyType QuadratureKeyType;

      typedef std::pair< MapperType, MapperType > MapperPairType;

    public:
      template <class QuadratureImpl>
      static const MapperPairType& getMapper(const QuadratureImpl& quadImpl,
                                             GeometryType elementGeometry,
                                             int faceIndex,
                                             int faceTwist)
      {
        // get quadrature implementation
        const QuadratureType& quad = quadImpl.ipList();

        // create key
        const QuadratureKeyType key (elementGeometry, quad.id() );

        MapperContainerType& mappers_ = mappers();
        MapperIteratorType it = mappers_.find( key );

        if( it == mappers_.end() )
        {
          std::integral_constant< bool, hasTwists > i2t;
          it = CacheProvider<GridPart, 1>::createMapper( quad, elementGeometry, i2t );
        }

        return it->second.getMapper(faceIndex, faceTwist);
      }

    private:
      typedef CacheStorage< ct, dim-codim, hasTwists>  CacheStorageType;
      typedef typename Traits::MapperVectorType MapperVectorType;

      typedef std::map<const QuadratureKeyType, CacheStorageType> MapperContainerType;
      typedef typename MapperContainerType::iterator MapperIteratorType;

    private:
      static MapperIteratorType
      createMapper ( const QuadratureType &quad, GeometryType elementGeometry, std::integral_constant< bool, true > );

      static MapperIteratorType
      createMapper ( const QuadratureType &quad, GeometryType elementGeometry, std::integral_constant< bool, false > );

      // mapper container holding mappings for different quads
      MapperContainerType mappers_;

    private:
      static MapperContainerType& mappers()
      {
        return instance().mappers_;
      }

      static ThisType& instance()
      {
        return Singleton< ThisType >::instance();
      }
    };

  } // namespace Fem

} // namespace Dune

#include "cacheprovider.cc"

#endif // #ifndef DUNE_FEM_CACHEPROVIDER_HH
