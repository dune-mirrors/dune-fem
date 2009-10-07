#ifndef DUNE_FEM_GENERICDOFMAPPER_HH
#define DUNE_FEM_GENERICDOFMAPPER_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referencetopologies.hh>

#include <dune/finiteelements/common/localcoefficients.hh>

#include <dune/fem/space/common/dofmapper.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridPart >
  class GenericDofMapper;

  template< class GridPart >
  class GenericDofMapIterator;



  // GenericDofMapperTraits
  // ----------------------

  template< class GridPart >
  struct GenericDofMapperTraits
  {
    typedef GenericDofMapper< GridPart > DofMapperType;

    typedef typename GridPart::GridType::template Codim< 0 >::Entity EntityType;

    typedef GenericDofMapIterator< GridPart > DofMapIteratorType;
  };



  // GenericDofMapper
  // ----------------

  template< class GridPart >
  class GenericDofMapper
  : public DofMapper< GenericDofMapperTraits< GridPart > >
  {
    typedef GenericDofMapper< GridPart > ThisType;
    typedef DofMapper< GenericDofMapperTraits< GridPart > > BaseType;

    friend class GenericDofMapIterator< GridPart >;

  public:
    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::DofMapIteratorType DofMapIteratorType;

    typedef GridPart GridPartType;

    struct SubEntityInfo
    {
      unsigned int codim;
      unsigned int subEntity;
      unsigned int topologyId;
      unsigned int numDofs;
      unsigned int offset;
    };

    struct MapInfo
    {
      typedef typename std::vector< SubEntityInfo >::const_iterator Iterator;

      Iterator begin() const 
      {
        return subEntityInfo.begin();
      }

      Iterator end() const 
      {
        return subEntityInfo.end();
      }

      unsigned int localDofPermutation( unsigned int i ) const
      {
        assert( i < localDof.size() );
        return localDof[ i ];
      }

      unsigned int size () const
      {
       return numDofs;
      }

      unsigned int numDofs;
      std::vector< unsigned int > localDof;
      std::vector< SubEntityInfo > subEntityInfo;
    };

  private:
    typedef typename GridPartType::IndexSetType IndexSetType;

    struct IndexInfo
    {
      unsigned int offset;
      unsigned int size;
      bool required;
      bool supported;

      IndexInfo () : supported( false ) {}
    };

    template< int topologyId >
    struct Build;

  public:
    static const int dimension = GridPartType::GridType::dimension;
    static const unsigned int numTopologies = (1 << dimension);

    template< class LocalCoefficientsProvider >
    GenericDofMapper ( const GridPartType &gridPart,
                       const LocalCoefficientsProvider &localCoefficientsProvider )
    : indexSet_( gridPart.indexSet() )
    {
      ForLoop< Build, 0, numTopologies-1 >::apply( localCoefficientsProvider, *this );
      update();
    }

    const IndexSetType &indexSet () const
    {
      return indexSet_;
    }

    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      typename DofMapIteratorType::Begin begin;
      return DofMapIteratorType( *this, entity, begin );
    }

    DofMapIteratorType end ( const EntityType &entity ) const
    {
      typename DofMapIteratorType::End end;
      return DofMapIteratorType( *this, entity, end );
    }

    template< class Entity >
    void map ( const Entity &entity, std::vector< unsigned int > &indices ) const;

    template< class Entity >
    int mapToGlobal( const Entity &entity, const int localDof ) const
    {
      static std::vector< unsigned int > indices;
      map( entity, indices );
      return indices[ localDof ];
    }

    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::mapEntityDofToGlobal not implemented, yet." );
    }

    unsigned int size ( const unsigned int topologyId ) const
    {
      assert( topologyId < numTopologies );
      return mapInfo_[ topologyId ].numDofs;
    }

    template< class Entity >
    unsigned int size ( const Entity &entity ) const
    {
      return size( GenericGeometry::topologyId( entity.type() ) );
    }

    unsigned int size () const
    {
      return size_;
    }

    bool topologyRequired ( const unsigned int topologyId ) const
    {
      const std::vector< IndexInfo > &indexInfo = indexInfo_[ 0 ];
      const unsigned int index                  = topologyId >> 1;
      assert( index < indexInfo.size() );
      return indexInfo[ index ].required;
    }

    int maxNumDofs () const
    {
      return maxNumDofs_;
    }

    template< class Entity >
    int numDofs ( const Entity &entity ) const
    {
      return mapInfo( entity ).numDofs;
    }

    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::numEntityDofs not implemented, yet." );
    }

    void update ();

    void update ( const bool overSizeMemory )
    {
      update();
    }

    template< class Entity >
    const MapInfo &mapInfo ( const Entity &entity ) const
    {
      return mapInfo_[ GenericGeometry::topologyId( entity.type() ) ];
    }

    const int numBlocks() const
    {
      return 1;
    }

    int offSet ( const int block ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::offSet not implemented, yet." );
    }

    int oldOffSet ( const int block ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::oldOffSet not implemented, yet." );
    }

    int numberOfHoles ( const int block ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::numberOfHoles not implemented, yet." );
    }

    int oldIndex ( const int hole, const int block) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::oldIndex not implemented, yet." );
    }

    int newIndex ( const int hole, const int block ) const
    {
      DUNE_THROW( NotImplemented, "GenericDofMapper::newIndex not implemented, yet." );
    }

    bool consecutive () const
    {
      return false;
    }

  private: 
    template< class Topology, class LocalCoefficients >
    void build ( const LocalCoefficients &localCoefficients );

    const IndexSetType &indexSet_;
    MapInfo mapInfo_[ numTopologies ];
    std::vector< IndexInfo > indexInfo_[ dimension+1 ];
    unsigned int maxNumDofs_;
    unsigned int size_;
  };


  template< class GridPart >
  template< class Entity >
  inline void GenericDofMapper< GridPart >
    ::map ( const Entity &entity, std::vector< unsigned int > &indices ) const
  {
    typedef typename std::vector< SubEntityInfo >::const_iterator Iterator;

    const MapInfo &mapInfo = mapInfo_[ GenericGeometry::topologyId( entity.type() ) ];

    indices.resize( mapInfo.numDofs );
    const unsigned int *localDof = &(mapInfo.localDof[ 0 ]);

    const Iterator end = mapInfo.subEntityInfo.end();
    for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
    {
      const unsigned int index  = indexSet().subIndex( entity, it->subEntity, it->codim );
      const unsigned int offset = it->offset + index*it->numDofs;
      for( unsigned int j = 0; j < it->numDofs; ++j )
        indices[ *(localDof++) ] = offset + j;
    }
  }


  template< class GridPart >
  inline void GenericDofMapper< GridPart >::update ()
  {
    maxNumDofs_ = 0;
    size_ = 0;
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int indexInfoSize = indexInfo_[ codim ].size();
      for( unsigned int i = 0; i < indexInfoSize; ++i )
      {
        IndexInfo &info = indexInfo_[ codim ][ i ];
        const unsigned int topologyId = i << 1;
        const unsigned int dim = dimension - codim;

        unsigned int idxSize = 0;
        if( GenericGeometry::hasGeometryType( topologyId, dim ) )
        {
          const GeometryType type = GenericGeometry::geometryType( topologyId, dim );
          idxSize = indexSet().size( type );
        }

        info.required = (idxSize > 0);
        assert( info.supported || !(info.required) );

        info.offset = size_;
        size_ += idxSize * info.size;
      }
    }

    for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      MapInfo &mapInfo = mapInfo_[ topologyId ];

      maxNumDofs_ = std::max( maxNumDofs_, mapInfo.numDofs );

      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
        it->offset = indexInfo_[ it->codim ][ it->topologyId >> 1 ].offset;
    }
  }


  template< class GridPart >
  template< class Topology, class LocalCoefficients >
  inline void GenericDofMapper< GridPart >
    ::build ( const LocalCoefficients &localCoefficients )
  {
    MapInfo &mapInfo = mapInfo_[ Topology::id ];

    mapInfo.numDofs = localCoefficients.size();
    mapInfo.localDof.resize( mapInfo.numDofs );

    const GenericGeometry::ReferenceTopology< dimension > &refTopology
      = GenericGeometry::ReferenceTopologies< dimension >::get( Topology::id );

    GenericGeometry::SubTopologyMapper< Topology > mapper;
    std::vector< unsigned int > counts( mapper.size(), (unsigned int)0 );

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      const LocalKey &key = localCoefficients.localKey( i );
      ++counts[ mapper( key.codim(), key.subEntity() ) ];
    }

    for( int codim = 0; codim <= dimension; ++codim )
    {
      const int subdimension = dimension-codim;
      // The last bit of the topology id is insignificat, hence store only
      // (1 << subdimension-1) many indexInfos.
      indexInfo_[ codim ].resize( subdimension > 0 ? 1 << (subdimension-1) : 1 );

      const unsigned int codimSize = refTopology.size( codim );
      for( unsigned int subEntity = 0; subEntity < codimSize; ++subEntity )
      {
        const unsigned int topologyId = refTopology.topologyId( codim, subEntity );
        IndexInfo &indexInfo          = indexInfo_[ codim ][ topologyId >> 1 ];

        const unsigned int count = counts[ mapper( codim, subEntity ) ];
        if( indexInfo.supported )
        {
          if( indexInfo.size != count )
            DUNE_THROW( InvalidStateException, "Inconsistent LocalCoefficients." );
        }
        else
          indexInfo.size = count;
        indexInfo.supported = true;

        if( count == 0 )
          continue;

        SubEntityInfo subEntityInfo;
        subEntityInfo.codim      = codim;
        subEntityInfo.subEntity  = subEntity;
        subEntityInfo.topologyId = topologyId;
        subEntityInfo.numDofs    = count;
        mapInfo.subEntityInfo.push_back( subEntityInfo );
      }
    }

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      const LocalKey &key = localCoefficients.localKey( i );

      unsigned int *localDof = &(mapInfo.localDof[ 0 ]);
      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); true; ++it )
      {
        if( it == end )
        {
          std::cerr << "Error: (subEntity = " << key.subEntity()
                    << ", codim = " << key.codim()
                    << ") not found in subEntityInfo" << std::endl;
          std::cerr << "SubEntityInfo contains:" << std::endl;
          for( it = mapInfo.subEntityInfo.begin(); it != end; ++it )
          {
            std::cerr << "  (subEntity = " << it->subEntity
                      << ", codim = " << it->codim << ")" << std::endl;
          }
          abort();
        }

        if( (it->codim == key.codim()) && (it->subEntity == key.subEntity()) )
        {
          *(localDof + key.index()) = i;
          break;
        }
        localDof += it->numDofs;
      }
    }
  }


  template< class GridPart >
  template< int topologyId >
  struct GenericDofMapper< GridPart >::Build
  {
    typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;

    template< class LocalCoefficientsProvider >
    static void
    apply ( const LocalCoefficientsProvider &localCoefficientsProvider,
            ThisType &dofMapper )
    {
      const unsigned int size = localCoefficientsProvider.template size< Topology >();
      switch( size )
      {
      case 0:
        break;

      case 1:
        dofMapper.template build< Topology >( localCoefficientsProvider.template localCoefficients< Topology >( 0 ) );
        break;

      default:
        DUNE_THROW( InvalidStateException, "GenericDofMapper supports at most one set of LocalCoefficients per topology." );
      }
    }
  };



  // GenericDofMapIterator
  // ---------------------
  
  template< class GridPart >
  class GenericDofMapIterator
  {
    typedef GenericDofMapIterator< GridPart > ThisType;

  public:
    typedef GenericDofMapper< GridPart > DofMapperType;

    typedef typename DofMapperType::EntityType EntityType;
    typedef typename DofMapperType::IndexSetType IndexSetType;

    typedef Int2Type< 0 > Begin;
    typedef Int2Type< 1 > End;

  private:
    typedef typename DofMapperType::SubEntityInfo SubEntityInfo;
    typedef typename DofMapperType::MapInfo MapInfo;
    typedef typename std::vector< SubEntityInfo >::const_iterator SubEntityInfoIterator;

  public:
    GenericDofMapIterator ( const DofMapperType &dofMapper,
                            const EntityType &entity,
                            const Begin &begin )
    : indexSet_( dofMapper.indexSet() ),
      entity_( entity ),
      mapInfo_( dofMapper.mapInfo_[ GenericGeometry::topologyId( entity.type() ) ] ),
      subEntityInfoIt_( mapInfo_.subEntityInfo.begin() ),
      localDof_( &(mapInfo_.localDof[ 0 ]) )
    {
      initSubEntity ( subEntityInfoIt_ );
    }

    GenericDofMapIterator ( const DofMapperType &dofMapper,
                            const EntityType &entity,
                            const End &end )
    : indexSet_( dofMapper.indexSet() ),
      entity_( entity ),
      mapInfo_( dofMapper.mapInfo_[ GenericGeometry::topologyId( entity.type() ) ] ),
      subEntityInfoIt_( mapInfo_.subEntityInfo.end() ),
      localDof_( &(mapInfo_.localDof[ mapInfo_.numDofs ]) )
    {}

    ThisType &operator++ ()
    {
      ++localDof_;
      ++dof_;
      if( dof_ >= subEntityInfoIt_->numDofs )
        initSubEntity( ++subEntityInfoIt_ );
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return (localDof_ == other.localDof_);
    }

    bool operator!= ( const ThisType &other ) const
    {
      return (localDof_ != other.localDof_);
    }

    int local () const
    {
      return *localDof_;
    }

    int global () const
    {
      return offset_ + dof_;
    }

  private:
    void initSubEntity ( const SubEntityInfoIterator &it )
    {
      if( it != mapInfo_.subEntityInfo.end() )
      {
        assert( it->numDofs > 0 );
        dof_ = 0;
        const unsigned int index = indexSet_.subIndex( entity_, it->subEntity, it->codim );
        offset_ = it->offset + index * it->numDofs;
      }
    }

  protected:
    const IndexSetType &indexSet_;
    const EntityType &entity_;
    const MapInfo &mapInfo_;
    SubEntityInfoIterator subEntityInfoIt_;
    const unsigned int *localDof_;
    int offset_;
    unsigned int dof_;
  };

}

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_GENERICDOFMAPPER_HH

/* vim: set sw=2 et: */
