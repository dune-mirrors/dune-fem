#ifndef DUNE_FEM_GENERICDOFMAPPER_HH
#define DUNE_FEM_GENERICDOFMAPPER_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referencetopologies.hh>

#include <dune/localfunctions/common/localkey.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridPart, class LocalCoefficientsMap >
  class GenericDofMapper;

  template< class GridPart, class LocalCoefficientsMap >
  class GenericDofMapIterator;



  // GenericDofMapperTraits
  // ----------------------

  template< class GridPart, class LocalCoefficientsMap >
  struct GenericDofMapperTraits
  {
    typedef GenericDofMapper< GridPart, LocalCoefficientsMap > DofMapperType;

    typedef typename GridPart::GridType::template Codim< 0 >::Entity EntityType;

    typedef GenericDofMapIterator< GridPart, LocalCoefficientsMap > DofMapIteratorType;
  };



  // GenericDofMapper
  // ----------------

  template< class GridPart, class LocalCoefficientsMap >
  class GenericDofMapper
  : public DofMapper< GenericDofMapperTraits< GridPart, LocalCoefficientsMap > >
  {
    typedef GenericDofMapper< GridPart, LocalCoefficientsMap > ThisType;
    typedef DofMapper< GenericDofMapperTraits< GridPart, LocalCoefficientsMap > > BaseType;

    friend class GenericDofMapIterator< GridPart, LocalCoefficientsMap >;

  public:
    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::DofMapIteratorType DofMapIteratorType;

    typedef GridPart GridPartType;
    typedef LocalCoefficientsMap LocalCoefficientsMapType;

    typedef typename LocalCoefficientsMapType::LocalCoefficientsType LocalCoefficientsType;

    struct SubEntityInfo
    {
      unsigned int codim;
      unsigned int subEntity;
      unsigned int topologyId;
      unsigned int blockIdx;
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
    typedef DofManager< typename GridPartType::GridType > DofManagerType;

    struct Block
    {
      unsigned int codim;
      unsigned int topologyId;
      unsigned int numDofs;
      unsigned int offset;
      unsigned int oldOffset;

      Block ( const unsigned int cd, const unsigned int tid, const unsigned int nDofs )
      : codim( cd ), topologyId( tid ), numDofs( nDofs )
      {}
    };

    template< int topologyId >
    struct Build;

  public:
    static const int dimension = GridPartType::GridType::dimension;
    static const unsigned int numTopologies = (1 << dimension);

    GenericDofMapper ( const GridPartType &gridPart,
                       const LocalCoefficientsMapType &localCoefficientsMap );

    ~GenericDofMapper ()
    {
      dofManager_.removeIndexSet( *this );
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
    const MapInfo &mapInfo ( const Entity &entity ) const;

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
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const;

    unsigned int size () const
    {
      return size_;
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
      const int codim = Entity::codimension;
      const unsigned int topologyId = entity.type().id(); // GenericGeometry::topologyId( entity.type() );
      const int blockIndex = blockIndex_[ codim ][ topologyId >> 1 ];
      return (blockIndex >= 0 ? blocks_[ blockIndex ].numDofs : 0);
    }

    bool contains ( const int codim ) const;

    const int numBlocks () const
    {
      return blocks_.size();
    }

    int offSet ( const int blockIdx ) const
    {
      assert( (blockIdx >= 0) && (blockIdx < numBlocks()) );
      return blocks_[ blockIdx ].offset;
    }

    int oldOffSet ( const int blockIdx ) const
    {
      assert( (blockIdx >= 0) && (blockIdx < numBlocks()) );
      return blocks_[ blockIdx ].oldOffset;
    }

    int numberOfHoles ( const int blockIdx ) const;
    int oldIndex ( const int hole, const int blockIdx ) const;
    int newIndex ( const int hole, const int blockIdx ) const;

    bool consecutive () const
    {
      return true;
    }

    bool fixedDataSize( const int codim ) const;


    // Adaptation Methods (as for Index Sets)

    template< class Entity >
    void insertEntity ( const Entity &entity )
    {
      update();
    }

    template< class Entity >
    void removeEntity ( const Entity &entity )
    {}

    void resize ()
    {
      update();
    }

    bool compress ()
    {
      update();
      return true;
    }

    void read_xdr ( const char *filename, int timestep )
    {
      update();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    GenericDofMapper ( const ThisType & );
    ThisType &operator= ( const ThisType & );

    void update ();

    template< class Topology >
    void build ( const LocalCoefficientsType &localCoefficients,
                 MapInfo &mapInfo );

    template< class Topology >
    void build ();

    DofManagerType &dofManager_;
    const IndexSetType &indexSet_;
    const LocalCoefficientsMapType &localCoefficientsMap_;
    std::vector< MapInfo > mapInfo_[ numTopologies ];
    std::vector< Block > blocks_;
    unsigned int maxNumDofs_;
    unsigned int size_;
    std::vector< int > blockIndex_[ dimension+1 ];
  };


  template< class GridPart, class LocalCoefficientsMap >
  inline GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::GenericDofMapper ( const GridPartType &gridPart,
                         const LocalCoefficientsMapType &localCoefficientsMap )
  : dofManager_( DofManagerType::instance( gridPart.grid() ) ),
    indexSet_( gridPart.indexSet() ),
    localCoefficientsMap_( localCoefficientsMap ),
    maxNumDofs_( 0 )
  {
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const int subdimension = dimension-codim;
      // The last bit of the topology id is insignificat, hence store only
      // (1 << subdimension-1) many indexInfos.
      blockIndex_[ codim ].resize( subdimension > 0 ? 1 << (subdimension-1) : 1, -1 );
    }
    ForLoop< Build, 0, numTopologies-1 >::apply( *this );
    update();
    dofManager_.addIndexSet( *this );
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline bool
  GenericDofMapper< GridPart, LocalCoefficientsMap >::contains ( const int codim ) const
  {
    typedef typename std::vector< int >::const_iterator Iterator;
    assert( (codim >= 0) && (codim <= dimension) );

    bool contains = false;
    const Iterator end = blockIndex_[ codim ].end();
    for( Iterator it  = blockIndex_[ codim ].begin(); it != end; ++it )
      contains |= (*it >= 0);
    return contains;
  }


  template< class GridPart, class LocalCoefficientsMap >
  template< class Entity >
  const typename GenericDofMapper< GridPart, LocalCoefficientsMap >::MapInfo &
  GenericDofMapper< GridPart, LocalCoefficientsMap >::mapInfo ( const Entity &entity ) const
  {
    const unsigned int topologyId = entity.type().id(); // GenericGeometry::topologyId( entity.type() );
    const unsigned int i = localCoefficientsMap_( entity );
    assert( i <= mapInfo_[ topologyId ].size() );
    return mapInfo_[ topologyId ][ i ];
  }


  template< class GridPart, class LocalCoefficientsMap >
  template< class Entity >
  inline void
  GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::map ( const Entity &entity, std::vector< unsigned int > &indices ) const
  {
    typedef typename std::vector< SubEntityInfo >::const_iterator Iterator;

    const MapInfo &info = mapInfo( entity );

    indices.resize( info.numDofs );
    const unsigned int *localDof = &(info.localDof[ 0 ]);

    const Iterator end = info.subEntityInfo.end();
    for( Iterator it = info.subEntityInfo.begin(); it != end; ++it )
    {
      const unsigned int index  = indexSet().subIndex( entity, it->subEntity, it->codim );
      const unsigned int offset = it->offset + index*it->numDofs;
      for( unsigned int j = 0; j < it->numDofs; ++j )
        indices[ *(localDof++) ] = offset + j;
    }
  }


  template< class GridPart, class LocalCoefficientsMap >
  template< class Entity >
  inline int
  GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
  {
    const int codim = Entity::codimension;
    const unsigned int topologyId = entity.type().id(); // GenericGeometry::topologyId( entity.type() );
    const int blockIndex = blockIndex_[ codim ][ topologyId >> 1 ];
    assert( blockIndex >= 0 );

    const Block &block = blocks_[ blockIndex ];
    assert( (unsigned int)localDof < block.numDofs );
    return block.offset + block.numDofs*indexSet_.index( entity ) + localDof;
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline void
  GenericDofMapper< GridPart, LocalCoefficientsMap >::update ()
  {
    size_ = 0;
    const unsigned int numBlocks = blocks_.size();
    for( unsigned int i = 0; i < numBlocks; ++i )
    {
      Block &block = blocks_[ i ];

      unsigned int idxSize = 0;
      const GeometryType type( block.topologyId, dimension - block.codim);
      if (!type.isNone())
        idxSize = indexSet().size( type );
      /*
      if( GenericGeometry::hasGeometryType( block.topologyId, dimension - block.codim ) )
      {
        const GeometryType type = GenericGeometry::geometryType( block.topologyId, dimension - block.codim );
        idxSize = indexSet().size( type );
      }
      */
      block.oldOffset = block.offset;
      block.offset = size_;
      size_ += idxSize * block.numDofs;
    }

    for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
    {
      typedef typename std::vector< MapInfo >::iterator Iterator;
      const Iterator mend = mapInfo_[ topologyId ].end();
      for( Iterator mit = mapInfo_[ topologyId ].begin(); mit != mend; ++mit )
      {
        typedef typename std::vector< SubEntityInfo >::iterator Iterator;
        const Iterator send = mit->subEntityInfo.end();
        for( Iterator sit = mit->subEntityInfo.begin(); sit != send; ++sit )
          sit->offset = blocks_[ sit->blockIdx ].offset;
      }
    }
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline int
  GenericDofMapper< GridPart, LocalCoefficientsMap >::numberOfHoles ( const int blockIdx ) const
  {
    assert( (blockIdx >= 0) && (blockIdx < numBlocks()) );
    const Block &block = blocks_[ blockIdx ];
    return block.numDofs * indexSet().numberOfHoles( block.codim );
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline int
  GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::oldIndex ( const int hole, const int blockIdx ) const
  {
    assert( (hole >= 0) && (hole < numberOfHoles( blockIdx )) );
    const Block &block = blocks_[ blockIdx ];
    const int numDofs = block.numDofs;
    return block.offset + numDofs * indexSet().oldIndex( hole / numDofs, block.codim ) + (hole % numDofs);
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline int
  GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::newIndex ( const int hole, const int blockIdx ) const
  {
    assert( (hole >= 0) && (hole < numberOfHoles( blockIdx )) );
    const Block &block = blocks_[ blockIdx ];
    const int numDofs = block.numDofs;
    return block.offset + numDofs * indexSet().newIndex( hole / numDofs, block.codim ) + (hole % numDofs);
  }


  template< class GridPart, class LocalCoefficientsMap >
  inline bool
  GenericDofMapper< GridPart, LocalCoefficientsMap >::fixedDataSize( const int codim ) const
  {
    typedef typename std::vector< int >::const_iterator Iterator;
    assert( (codim >= 0) && (codim <= dimension) );

    Iterator begin = blockIndex_[ codim ].begin();
    const Iterator end = blockIndex_[ codim ].end();
    if( begin == end )
      return true;

    unsigned int numDofs = (*begin >= 0 ? blocks_[ *begin ].numDofs : 0);
    bool fixedSize = true;
    for( Iterator it = begin++; it != end; ++it )
      fixedSize &= (numDofs == (*it >= 0 ? blocks_[ *it ].numDofs : 0));
    return fixedSize;
  }


  template< class GridPart, class LocalCoefficientsMap >
  template< class Topology >
  inline void
  GenericDofMapper< GridPart, LocalCoefficientsMap >
    ::build ( const LocalCoefficientsType &localCoefficients,
              MapInfo &mapInfo )
  {
    const GenericGeometry::ReferenceTopology< dimension > &refTopology
      = GenericGeometry::ReferenceTopologies< dimension >::get( Topology::id );

    mapInfo.numDofs = localCoefficients.size();
    mapInfo.localDof.resize( mapInfo.numDofs );
    maxNumDofs_ = std::max( maxNumDofs_, mapInfo.numDofs );

    // count the number of DoFs on each subentity
    GenericGeometry::SubTopologyMapper< Topology > mapper;
    std::vector< unsigned int > counts( mapper.size(), (unsigned int)0 );
    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      const LocalKey &key = localCoefficients.localKey( i );
      ++counts[ mapper( key.codim(), key.subEntity() ) ];
    }

    // build subentity information
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int codimSize = refTopology.size( codim );
      for( unsigned int subEntity = 0; subEntity < codimSize; ++subEntity )
      {
        const unsigned int topologyId = refTopology.topologyId( codim, subEntity );

        int &blockIdx = blockIndex_[ codim ][ topologyId >> 1 ];
        const unsigned int numDofs = counts[ mapper( codim, subEntity ) ];
        if( blockIdx == -1 )
        {
          if( numDofs > 0 )
          {
            blockIdx = blocks_.size();
            blocks_.push_back( Block( codim, topologyId, numDofs ) );
          }
          else
            blockIdx = -2;
        }
        else if( numDofs != (blockIdx >= 0 ? blocks_[ blockIdx ].numDofs : 0) )
          DUNE_THROW( InvalidStateException, "Inconsistent LocalCoefficients." );

        if( numDofs > 0 )
        {
          SubEntityInfo subEntityInfo;
          subEntityInfo.codim      = codim;
          subEntityInfo.subEntity  = subEntity;
          subEntityInfo.topologyId = topologyId;
          subEntityInfo.blockIdx   = blockIdx;
          subEntityInfo.numDofs    = numDofs;
          mapInfo.subEntityInfo.push_back( subEntityInfo );
        }
      }
    }

    // build permutation of local dofs
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


  template< class GridPart, class LocalCoefficientsMap >
  template< class Topology >
  inline void
  GenericDofMapper< GridPart, LocalCoefficientsMap >::build ()
  {
    const unsigned int size = localCoefficientsMap_.template size< Topology >();
    mapInfo_[ Topology::id ].resize( size );
    for( unsigned int i = 0; i < size; ++i )
    {
      MapInfo &mapInfo = mapInfo_[ Topology::id ][ i ];
      const LocalCoefficientsType &localCoefficients
        = localCoefficientsMap_.template localCoefficients< Topology >( i );
      build< Topology >( localCoefficients, mapInfo );
    }
  }



  template< class GridPart, class LocalCoefficientsMap >
  template< int topologyId >
  struct GenericDofMapper< GridPart, LocalCoefficientsMap >::Build
  {
    typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;

    static void apply ( ThisType &dofMapper )
    {
      dofMapper.template build< Topology >();
    }
  };



  // GenericDofMapIterator
  // ---------------------
  
  template< class GridPart, class LocalCoefficientsMap >
  class GenericDofMapIterator
  {
    typedef GenericDofMapIterator< GridPart, LocalCoefficientsMap > ThisType;

  public:
    typedef GenericDofMapper< GridPart, LocalCoefficientsMap > DofMapperType;

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
      mapInfo_( dofMapper.mapInfo( entity ) ),
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
      mapInfo_( dofMapper.mapInfo( entity ) ),
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
