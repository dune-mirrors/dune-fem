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

    template< class LocalCoefficientsProvider >
    GenericDofMapper ( const GridPartType &gridPart,
                       const LocalCoefficientsProvider &localCoefficientsProvider );

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
      const int codim = Entity::codimension;
      const unsigned int topologyId = GenericGeometry::topologyId( entity.type() );
      const int blockIndex = blockIndex_[ codim ][ topologyId >> 1 ];
      assert( blockIndex >= 0 );

      const Block &block = blocks_[ blockIndex ];
      assert( (unsigned int)localDof < block.numDofs );
      return block.offset + block.numDofs*indexSet_.index( entity ) + localDof;
    }

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
      const unsigned int topologyId = GenericGeometry::topologyId( entity.type() );
      const int blockIndex = blockIndex_[ codim ][ topologyId >> 1 ];
      return (blockIndex >= 0 ? blocks_[ blockIndex ].numDofs : 0);
    }

    bool contains ( const int codim ) const
    {
      typedef typename std::vector< int >::const_iterator Iterator;
      assert( (codim >= 0) && (codim <= dimension) );

      bool contains = false;
      const Iterator end = blockIndex_[ codim ].end();
      for( Iterator it  = blockIndex_[ codim ].begin(); it != end; ++it )
        contains |= (*it >= 0);
      return contains;
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

    int numberOfHoles ( const int blockIdx ) const
    {
      assert( (blockIdx >= 0) && (blockIdx < numBlocks()) );
      const Block &block = blocks_[ blockIdx ];
      return block.numDofs * indexSet().numberOfHoles( block.codim );
    }

    int oldIndex ( const int hole, const int blockIdx ) const
    {
      assert( (hole >= 0) && (hole < numberOfHoles( blockIdx )) );
      const Block &block = blocks_[ blockIdx ];
      const int numDofs = block.numDofs;
      return block.offset + numDofs * indexSet().oldIndex( hole / numDofs, block.codim ) + (hole % numDofs);
    }

    int newIndex ( const int hole, const int blockIdx ) const
    {
      assert( (hole >= 0) && (hole < numberOfHoles( blockIdx )) );
      const Block &block = blocks_[ blockIdx ];
      const int numDofs = block.numDofs;
      return block.offset + numDofs * indexSet().newIndex( hole / numDofs, block.codim ) + (hole % numDofs);
    }

    bool consecutive () const
    {
      return true;
    }

    bool fixedDataSize( const int codim ) const
    {
      typedef typename std::vector< int >::const_iterator Iterator;
      assert( (codim >= 0) && (codim <= dimension) );

      const Iterator begin = blockIndex_[ codim ].begin();
      const Iterator end = blockIndex_[ codim ].end();
      if( begin == end )
        return true;

      unsigned int numDofs = (*begin >= 0 ? blocks_[ *begin ].numDofs : 0);
      bool fixedSize = true;
      for( Iterator it = begin++; it != end; ++it )
        fixedSize &= (numDofs == (*it >= 0 ? blocks_[ *it ].numDofs : 0));
      return fixedSize;
    }

  private: 
    template< class Topology, class LocalCoefficients >
    void build ( const LocalCoefficients &localCoefficients,
                 MapInfo &mapInfo );

    const IndexSetType &indexSet_;
    MapInfo mapInfo_[ numTopologies ];
    std::vector< Block > blocks_;
    unsigned int maxNumDofs_;
    unsigned int size_;
    std::vector< int > blockIndex_[ dimension+1 ];
  };


  template< class GridPart >
  template< class LocalCoefficientsProvider >
  inline GenericDofMapper< GridPart >
    ::GenericDofMapper ( const GridPartType &gridPart,
                         const LocalCoefficientsProvider &localCoefficientsProvider )
  : indexSet_( gridPart.indexSet() ),
    maxNumDofs_( 0 )
  {
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const int subdimension = dimension-codim;
      // The last bit of the topology id is insignificat, hence store only
      // (1 << subdimension-1) many indexInfos.
      blockIndex_[ codim ].resize( subdimension > 0 ? 1 << (subdimension-1) : 1, -1 );
    }
    ForLoop< Build, 0, numTopologies-1 >::apply( localCoefficientsProvider, *this );
    update();
  }


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
    size_ = 0;
    const unsigned int numBlocks = blocks_.size();
    for( unsigned int i = 0; i < numBlocks; ++i )
    {
      Block &block = blocks_[ i ];

      unsigned int idxSize = 0;
      if( GenericGeometry::hasGeometryType( block.topologyId, dimension - block.codim ) )
      {
        const GeometryType type = GenericGeometry::geometryType( block.topologyId, dimension - block.codim );
        idxSize = indexSet().size( type );
      }

      block.oldOffset = block.offset;
      block.offset = size_;
      size_ += idxSize * block.numDofs;
    }

    for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      MapInfo &mapInfo = mapInfo_[ topologyId ];

      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
        it->offset = blocks_[ it->blockIdx ].offset;
    }
  }


  template< class GridPart >
  template< class Topology, class LocalCoefficients >
  inline void GenericDofMapper< GridPart >
    ::build ( const LocalCoefficients &localCoefficients,
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
      if( size > 1 )
        DUNE_THROW( InvalidStateException, "GenericDofMapper supports at most one set of LocalCoefficients per topology." );
      for( unsigned int i = 0; i < size; ++i )
      {
        MapInfo &mapInfo = dofMapper.mapInfo_[ topologyId ];
        dofMapper.template build< Topology >( localCoefficientsProvider.template localCoefficients< Topology >( 0 ), mapInfo );
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
