#ifndef DUNE_FEM_PERIODICINDEXSET_HH
#define DUNE_FEM_PERIODICINDEXSET_HH

//- system includes
#include <vector>

//- dune inclues
#include <dune/grid/common/grid.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>

namespace Dune
{

  //! This index set supports only codimensions 0 and dimension!
  template< class Grid >
  class PeriodicLeafIndexSet
  : public DuneGridIndexSetAdapter
    < Grid, PeriodicLeafIndexSet< Grid >, DefaultLeafIteratorTypes< Grid > >
  {
    typedef DefaultLeafIteratorTypes< Grid > IndexSetTypes;
    
    typedef PeriodicLeafIndexSet< Grid > ThisType;
    typedef DuneGridIndexSetAdapter< Grid, ThisType, IndexSetTypes > BaseType;

    friend class DuneGridIndexSetAdapter< Grid, ThisType, IndexSetTypes >;
    
  public:
    //! type of the grid
    typedef Grid GridType;

    //! type of indices
    typedef typename BaseType :: IndexType IndexType;

    //! coordinate type of the grid
    typedef typename GridType :: ctype ctype;
     
    //! dimension of the grid
    enum { dimension = GridType :: dimension };

  protected:
    typedef typename GridType :: Traits :: LeafIndexSet BaseIndexSetType;

    static const ctype tolerance = 1e-10;

#ifdef INDEXSET_HAS_ITERATORS
  public:
    template< int codim >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename IndexSetTypes :: template Codim< codim >
          :: template Partition< pitype > :: Iterator
          Iterator;
      };
    };
#endif

  protected:
    using BaseType :: grid_;
    const BaseIndexSetType &baseIndexSet_;
    
    IndexType size_;
    DynamicArray< int > index_;
    
    IndexType edgeSize_;
    DynamicArray< int > edgeIndex_;

  public:
    //! Create PeriodicIndexSet for a grid
    inline explicit PeriodicLeafIndexSet ( const GridType &grid );

    template< class Entity >
    inline IndexType index ( const Entity &entity ) const
    {
      return index< Entity :: codimension >( entity );
    }
    
    template< int codim >
    inline IndexType
    index ( const typename GridType :: template Codim< codim > :: Entity &entity ) const
    {
      return map< codim >( baseIndexSet_.template index< codim >( entity ) );
    }

    template< int codim >
    inline IndexType
    subIndex ( const typename GridType :: template Codim< 0 > ::  Entity &entity,
               int subNumber ) const
    {
      return map< codim >( baseIndexSet_.template subIndex< codim >( entity, subNumber ) );
    }

    //! Return true, if this set provides an index for an entity
    template< class Entity >
    inline bool contains ( const Entity &entity ) const;

    //! Return the size of the index set for a codimension
    inline IndexType size ( int codim ) const;

    //! Return the size of the index set for a codimension
    //! Marked for revision in grid/common/defaultindexsets.hh
    inline IndexType size ( GeometryType type ) const
    {
      if( geometryTypeValid( type ) )
        return size( dimension - type.dim() );
      return 0;
    }

    //! Deliver all geometry types used in this grid
    inline const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      return baseIndexSet_.geomTypes( codim );
    }

#ifdef INDEXSET_HAS_ITERATORS
    /** \brief Deliver iterator to first entity of given codimension and partition type.
     */
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: Iterator begin () const
    {
      return grid_.template leafbegin< codim, pitype >();
    }

    /** \brief Deliver iterator to last entity of given codimension and partition type.
     */
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: Iterator end () const
    {
      return grid_.template leafend< codim, pitype >();
    }
#endif

    inline std :: string name () const
    {
      return std :: string( "PeriodicLeafIndexSet" );
    }

  protected:
    template< class field, int dim >
    inline static field
    dist2 ( const FieldVector< field, dim > &u, const FieldVector< field, dim > &v );

    template< class field, int dim >
    inline static void
    truncate ( FieldVector< field, dim > &v );
    template< class field, int dim >
    inline static void
    truncate ( FieldVector< field, dim > &v, FieldVector< field, dim > &w );
   
    template< class Iterator >
    inline void calcVertexIndices ( const Iterator &begin, const Iterator &end );
    template< class Iterator >
    inline void calcEdgeIndices ( const Iterator &begin, const Iterator &end );
    
    bool geometryTypeValid( const GeometryType &type ) const;

    template< int codim >
    inline IndexType map ( const IndexType &baseIndex ) const;
  };


    
  template< class Grid >
  inline PeriodicLeafIndexSet< Grid >
    :: PeriodicLeafIndexSet ( const GridType &grid )
  : BaseType( grid ),
    baseIndexSet_( grid.leafIndexSet() )
  {
    calcVertexIndices( grid.template leafbegin< dimension, All_Partition >(),
                       grid.template leafend< dimension, All_Partition >() );
    if( dimension > 1 )
    {
      const int codim = dimension - 1;
      calcEdgeIndices( grid.template leafbegin< codim, All_Partition >(),
                       grid.template leafend< codim, All_Partition >() );
    }
  }

  template< class Grid >
  template< class Entity >
  inline bool
  PeriodicLeafIndexSet< Grid > :: contains ( const Entity &entity ) const
  {
    if( !baseIndexSet_.contains( entity ) )
      return false;

    switch( Entity :: codimension )
    {
    case 0:
      return true;

    case dimension - 1:
      return edgeIndex_[ baseIndexSet_.index( entity ) ] >= 0;

    case dimension:
      return index_[ baseIndexSet_.index( entity ) ] >= 0;

    default:
      return false;
    }
  }

  template< class Grid >
  inline typename PeriodicLeafIndexSet< Grid > :: IndexType
  PeriodicLeafIndexSet< Grid > :: size ( int codim ) const
  {
    switch( codim )
    {
    case 0:
      return baseIndexSet_.size( 0 );

    case dimension - 1:
      return edgeSize_;

    case dimension:
      return size_;
    default:
      return 0;
    }
  }

  template< class Grid >
  template< class field, int dim >
  inline field
  PeriodicLeafIndexSet< Grid >
    :: dist2 ( const FieldVector< field, dim > &u, const FieldVector< field, dim > &v )
  {
    field d2 = 0;
    for( int i = 0; i < dim; ++i )
    {
      const field d = u[ i ] - v[ i ];
      d2 += d * d;
    }
    return d2;
  }

  template< class Grid >
  template< class field, int dim >
  inline void
  PeriodicLeafIndexSet< Grid > :: truncate ( FieldVector< field, dim > &v )
  {
     for( int i = 0; i < dim; ++i )
       v[ i ] -= floor( v[ i ] );
  }

  template< class Grid >
  template< class field, int dim >
  inline void
  PeriodicLeafIndexSet< Grid >
    :: truncate ( FieldVector< field, dim > &v, FieldVector< field, dim > &w )
  {
    for( int i = 0; i < dim; ++i )
    {
      const field f = std :: min( floor( v[ i ] ), floor( w[ i ] ) );
      v[ i ] -= f;
      w[ i ] -= f;
    }
  }

  template< class Grid >
  template< class Iterator >
  inline void
  PeriodicLeafIndexSet< Grid >
    :: calcVertexIndices ( const Iterator &begin, const Iterator &end )
  {
    enum { codim = dimension };

    typedef typename GridType :: template Codim< codim > :: Geometry GeometryType;
    const unsigned int dimWorld = GeometryType :: coorddimension;
    typedef FieldVector< typename GeometryType :: ctype, dimWorld > VertexType;
    
    const int n = baseIndexSet_.size( codim );
    index_.resize( n );
    index_.assign( -1 );
    size_ = 0;
    
    DynamicArray< VertexType > vertices( n );
    
    for( Iterator it = begin; it != end; ++it )
    {
      int idx = baseIndexSet_.index( *it );
      VertexType &vertex = vertices[ idx ];
      
      const GeometryType &geometry = it->geometry();
      vertex = geometry[ 0 ];
      truncate( vertex );

      for( int i = 0; i < n; ++i )
      {
        if( index_[ i ] < 0 )
          continue;

        if( dist2( vertex, vertices[ i ] ) < tolerance*tolerance )
        {
          index_[ idx ] = index_[ i ];
          break;
        }
      }

      if( index_[ idx ] < 0 )
        index_[ idx ] = size_++;
    }
  }

  template< class Grid >
  template< class Iterator >
  inline void
  PeriodicLeafIndexSet< Grid >
    :: calcEdgeIndices ( const Iterator &begin, const Iterator &end )
  {
    enum { codim = dimension - 1 };

    typedef typename GridType :: template Codim< codim > :: Geometry GeometryType;
    typedef FieldVector< typename GeometryType :: ctype,
                         GeometryType :: coorddimension >
      VertexType;

    const ctype eps = tolerance * tolerance;
    
    const int n = baseIndexSet_.size( codim );
    edgeIndex_.resize( n );
    edgeIndex_.assign( -1 );
    edgeSize_ = 0;

    DynamicArray< VertexType > vertices0( n );
    DynamicArray< VertexType > vertices1( n );
    
    for( Iterator it = begin; it != end; ++it )
    {
      int idx = baseIndexSet_.index( *it );
      VertexType &vertex0 = vertices0[ idx ];
      VertexType &vertex1 = vertices1[ idx ];
      
      const GeometryType &geometry = it->geometry();
      vertex0 = geometry[ 0 ];
      vertex1 = geometry[ 1 ];
      truncate( vertex0, vertex1 );

      for( int i = 0; i < n; ++i )
      {
        if( edgeIndex_[ i ] < 0 )
          continue;

        bool same = (dist2( vertex0, vertices0[ i ] ) < eps)
                    && (dist2( vertex1, vertices1[ i ] ) < eps);
        same |= (dist2( vertex0, vertices1[ i ] ) < eps)
                && (dist2( vertex1, vertices0[ i ] ) < eps);
        if( same )
        {
          edgeIndex_[ idx ] = edgeIndex_[ i ];
          break;
        }
      }

      if( edgeIndex_[ idx ] < 0 )
        edgeIndex_[ idx ] = edgeSize_++;
    }
  }

  template< class Grid >
  inline bool
  PeriodicLeafIndexSet< Grid > :: geometryTypeValid( const GeometryType &type ) const
  {
    const int codim = dimension - type.dim();
    const std :: vector< GeometryType > &gT = geomTypes( codim );
    for( size_t i = 0; i < gT.size(); ++i )
    {
      if( gT[ i ] == type )
        return true;
    }
    return false;
  }

  template< class Grid >
  template< int codim >
  inline typename PeriodicLeafIndexSet< Grid > :: IndexType
  PeriodicLeafIndexSet< Grid > :: map ( const IndexType &baseIndex ) const
  {
    switch( codim )
    {
    case 0:
      return baseIndex;

    case dimension - 1:
      assert( edgeIndex_[ baseIndex ] >= 0 );
      return edgeIndex_[ baseIndex ];

    case dimension:
      assert( index_[ baseIndex ] >= 0 );
      return index_[ baseIndex ];

    default:
      DUNE_THROW( NotImplemented, "PeriodicLeafIndexSet supports only codimensions 0, (dimension-1) and dimension!" );
    }
    return 0;
  }

}

#endif
