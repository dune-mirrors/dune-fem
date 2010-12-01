#ifndef DUNE_CACHEPROVIDER_HH
#define DUNE_CACHEPROVIDER_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/common/misc.hh>
#include <dune/common/version.hh>

#include <dune/fem/misc/capabilities.hh>

//- Local includes
#include "pointmapper.hh"
#include "twistprovider.hh"
#include "pointprovider.hh"

namespace Dune {

  //! Storage class for mappers.
  template <class ct, int dim, bool unstructured>
  class CacheStorage;

  //! Specialisation for grids with twist (i.e. unstructured ones).
  template <class ct, int dim>
  class CacheStorage<ct, dim, true>
  {
  private:
    typedef CachingTraits<ct, dim> Traits;

  public:
    typedef typename Traits::MapperType MapperType;
    
  public:
    CacheStorage(int numFaces, int maxTwist) :
      mappers_(numFaces)
    {
      for (MapperIteratorType it = mappers_.begin(); 
           it != mappers_.end(); ++it) {
        it->resize(maxTwist + Traits::twistOffset_);
      }
    }

    CacheStorage(const CacheStorage& other) :
      mappers_(other.mappers_)
    {}

    void addMapper(const MapperType& faceMapper, 
                   const MapperType& twistMapper,
                   int faceIndex, int faceTwist)
    {
      assert(twistMapper.size() == faceMapper.size());

      MapperType& mapper = 
        mappers_[faceIndex][faceTwist + Traits::twistOffset_];
      mapper.resize(twistMapper.size());

      for (size_t i = 0; i < mapper.size(); ++i) {
        mapper[i] = faceMapper[twistMapper[i]];
      }

    }

    const MapperType& getMapper(int faceIndex, int faceTwist) const
    {
      assert( faceTwist + Traits::twistOffset_ >= 0 );
      return mappers_[faceIndex][faceTwist + Traits::twistOffset_];
    }

  private:
    typedef std::vector<std::vector<MapperType> > MapperContainerType;
    typedef typename MapperContainerType::iterator MapperIteratorType;

  private:
    MapperContainerType mappers_;
  };

  template <class ct, int dim>
  class CacheStorage<ct, dim, false>
  {
  private:
    typedef CachingTraits<ct, dim> Traits;

  public:
    typedef typename Traits::MapperType MapperType;

  public:
    explicit CacheStorage ( int numFaces )
    : mappers_( numFaces )
    {}

    CacheStorage ( const CacheStorage &other )
    : mappers_( other.mappers_ )
    {}

    void addMapper ( const MapperType &mapper, int faceIndex )
    {
      assert( (faceIndex >= 0) && (faceIndex < (int)mappers_.size()) );
      mappers_[ faceIndex ] = mapper;
    }

    const MapperType &getMapper ( int faceIndex, int faceTwist ) const
    {
#ifndef NDEBUG 
      if( faceIndex >= (int)mappers_.size() )
        std::cerr << "Error: faceIndex = " << faceIndex << " >= " << mappers_.size() << " = mappers_.size()" << std::endl;
#endif
      assert( (faceIndex >= 0) && (faceIndex < (int)mappers_.size()) );
      return mappers_[ faceIndex ];
    }

  private:
    typedef typename std::vector<MapperType> MapperContainerType;

  private:
    MapperContainerType mappers_;
  };


  // CacheProvider
  // -------------

  template< class GridImp, int codim >
  class CacheProvider;

  template <class GridImp>
  class CacheProvider<GridImp, 0>
  {
  private:
    enum { codim = 0 };
    enum { dim = GridImp::dimension };
    typedef typename GridImp::ctype ct;
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

  template <class GridImp>
  class CacheProvider<GridImp, 1>
  {
    enum { codim = 1 };
    enum { dim = GridImp::dimension };
    typedef typename GridImp::ctype ct;
    typedef CachingTraits<ct, dim-codim> Traits;

  public:
    typedef typename Traits::QuadratureType QuadratureType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::QuadratureKeyType QuadratureKeyType;

  public:
    template <class QuadratureImpl>
    static const MapperType& getMapper(const QuadratureImpl& quadImpl,
                                       GeometryType elementGeometry,
                                       int faceIndex,
                                       int faceTwist)
    {
      // get quadrature implementation 
      const QuadratureType& quad = quadImpl.ipList();
      
      // create key 
      const QuadratureKeyType key (elementGeometry, quad.id() );
      
      MapperIteratorType it = mappers_.find( key );
      
      if( it == mappers_.end() ) 
      {
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,1,0)
        integral_constant< bool, Capabilities::IsUnstructured< GridImp >::v > i2t;
#else
        Int2Type< Capabilities::IsUnstructured< GridImp >::v > i2t;
#endif
        it = CacheProvider<GridImp, 1>::createMapper( quad, elementGeometry, i2t );
      }

      return it->second.getMapper(faceIndex, faceTwist);
    }

  private:
    typedef CacheStorage<
      ct, dim-codim,
       Capabilities::IsUnstructured<GridImp>::v> 
        CacheStorageType; 
    typedef typename Traits::MapperVectorType MapperVectorType;
    typedef std::map<const QuadratureKeyType, CacheStorageType> MapperContainerType;
    typedef typename MapperContainerType::iterator MapperIteratorType;

  private:
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,1,0)
    static MapperIteratorType
    createMapper ( const QuadratureType &quad, GeometryType elementGeometry, integral_constant< bool, true > );

    static MapperIteratorType
    createMapper ( const QuadratureType &quad, GeometryType elementGeometry, integral_constant< bool, false > );
#else
    static MapperIteratorType
    createMapper ( const QuadratureType &quad, GeometryType elementGeometry, Int2Type< true > );

    static MapperIteratorType
    createMapper ( const QuadratureType &quad, GeometryType elementGeometry, Int2Type< false > );
#endif

  private:
    static MapperContainerType mappers_;
  };

}

#include "cacheprovider.cc"

#endif // #ifndef DUNE_CACHEPROVIDER_HH
