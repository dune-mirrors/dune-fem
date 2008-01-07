#ifndef DUNE_FEM_PERIODICINDEXSET_HH
#define DUNE_FEM_PERIODICINDEXSET_HH

//- system includes
#include <vector>

//- dune inclues
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/fem/storage/array.hh>

namespace Dune
{

  //! This index set supports only codimensions 0 and dimension!
  template< class GridImp >
  class PeriodicLeafIndexSet
  : public IndexSetDefaultImplementation
    < GridImp, PeriodicLeafIndexSet< GridImp >, DefaultLeafIteratorTypes< GridImp > >
  {
  public:
    //! type of the grid
    typedef GridImp GridType;
      
    //! coordinate type of the grid
    typedef typename GridType :: ctype ctype;
     
    //! dimension of the grid
    enum { dimension = GridType :: dimension };

  private:
    typedef DefaultLeafIteratorTypes< GridType > IndexSetTypes;

    typedef PeriodicLeafIndexSet< GridType > ThisType;
    typedef IndexSetDefaultImplementation< GridType, ThisType, IndexSetTypes >
      BaseType;

  protected:
    typedef typename GridType :: Traits :: LeafIndexSet BaseIndexSetType;

    static const ctype tolerance = 1e-10;

  public:
    template< int codim >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename IndexSetTypes :: template Codim< codim >
                                       :: template Partition< pitype >
                                       :: Iterator
          Iterator;
      };
    };

  protected:
    const GridType &grid_;
    const BaseIndexSetType &baseIndexSet_;
    
    int size_;
    DynamicArray< int > index_;
    
    int edgeSize_;
    DynamicArray< int > edgeIndex_;

  public:
    //! Create PeriodicIndexSet for a grid
    inline explicit PeriodicLeafIndexSet ( const GridType &grid )
    : grid_( grid ),
      baseIndexSet_( grid.leafIndexSet() )
    {
      calcVertexIndices();
      if( dimension > 1 )
        calcEdgeIndices();
    }

    //! Return periodic index of an entity
    template< class EntityType >
    int index ( const EntityType &entity ) const
    {
      return index< EntityType :: codimension >( entity );
    }

    template< int codim >
    int index ( const typename GridType :: template Codim< codim > :: Entity &entity ) const
    {
      const int baseIndex = baseIndexSet_.index( entity );
      
      if( codim == 0 )
        return baseIndex;

      if( codim == dimension )
      {
        assert( index_[ baseIndex ] >= 0 );
        return index_[ baseIndex ];
      }

      if( codim == dimension - 1 )
      {
        assert( edgeIndex_[ baseIndex ] >= 0 );
        return edgeIndex_[ baseIndex ];
      }
      
      DUNE_THROW( NotImplemented, "PeriodicLeafIndexSet supports only codimensions 0, (dimension-1) and dimension!" );
      return 0;
    }

    //! Given a codimension and a number, return the associated subentity of
    //! an entity of codimension 0.
    template< int codim >
    int subIndex ( const typename GridType :: template Codim< 0 > :: Entity &entity, int i ) const
    {
      const int baseIndex = baseIndexSet_.template subIndex< codim >( entity, i );
      
      if( codim == 0 )
        return baseIndex;
      
      if( codim == dimension )
      {
        assert( index_[ baseIndex ] >= 0 );
        return index_[ baseIndex ];
      }

      if( codim == dimension - 1 )
      {
        assert( edgeIndex_[ baseIndex ] >= 0 );
        return edgeIndex_[ baseIndex ];
      }

      DUNE_THROW( NotImplemented, "PeriodicLeafIndexSet supports only codimensions 0, (dimension-1) and dimension!" );
      return 0;
    }

    //! Return true, if this set provides an index for an entity
    template< class EntityType >
    bool contains ( const EntityType &entity ) const
    {
      enum { codim = EntityType :: codimension };

      if( codim == 0 )
        return true;

      if( codim == dimension )
        return index_[ baseIndexSet_.index( entity ) ] >= 0;
      
      if( codim == dimension - 1 )
        return edgeIndex_[ baseIndexSet_.index( entity ) ] >= 0;

      return false;
    }

    //! Return the size of the index set for a codimension
    int size ( int codim ) const
    {
      if( codim == 0)
        return baseIndexSet_.size( 0 );

      if( codim == dimension )
        return size_;

      if( codim == dimension - 1 )
        return edgeSize_;
 
      return 0;
    }

    //! Return the size of the index set for a codimension
    //! Marked for revision in grid/common/defaultindexsets.hh
    int size ( GeometryType type ) const
    {
      if( geometryTypeValid( type ) )
        return size( dimension - type.dim() );
      return 0;
    }

    //! Deliver all geometry types used in this grid
    const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      return baseIndexSet_.geomTypes( codim );
    }

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

    inline const std :: string name () const
    {
      return std :: string( "PeriodicLeafIndexSet" );
    }

  protected:
    // Note: This should be placed in FieldVector!
    //       Moreover, this could be implemented via template metaprogramming
    template< class ctype, int dim >
    ctype dist2 ( const FieldVector< ctype, dim > &u, const FieldVector< ctype, dim > &v )
    {
      ctype d2 = 0;
      for( int i = 0; i < dim; i++ ) {
          ctype d = u[ i ] - v[ i ];
          d2 += d * d;
      }
      return d2;
    }

    template< class ctype, int dim >
    void truncate ( FieldVector< ctype, dim > &v )
    {
       for( int i = 0; i < dim; ++i )
         v[ i ] -= floor( v[ i ] );
    }

    template< class ctype, int dim >
    void truncate ( FieldVector< ctype, dim > &v, FieldVector< ctype, dim > &w )
    {
      for( int i = 0; i < dim; ++i ) {
          ctype f = fmin( floor( v[ i ] ), floor( w[ i ] ) );
          v[ i ] -= f;
          w[ i ] -= f;
      }
    }
    
    void calcVertexIndices ()
    {
      enum { codim = dimension };

      typedef typename Codim< codim > :: template Partition< All_Partition > :: Iterator
        IteratorType;

      typedef typename GridType :: template Codim< codim > :: Geometry GeometryType;
      typedef FieldVector< typename GeometryType :: ctype, 
                           GeometryType :: coorddimension >
        VertexType;

      
      const int n = baseIndexSet_.size( codim );
      index_.resize( n );
      index_.assign( -1 );
      size_ = 0;
      
      DynamicArray< VertexType > vertices( n );
      
      const IteratorType eit = end< codim, All_Partition >();
      for( IteratorType it = begin< codim, All_Partition >(); it != eit; ++it )
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

    void calcEdgeIndices ()
    {
      enum { codim = dimension - 1 };

      typedef typename Codim< codim > :: template Partition< All_Partition > :: Iterator
        IteratorType;
      
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
      
      IteratorType eit = end< codim, All_Partition >();
      for( IteratorType it = begin< codim, All_Partition >(); it != eit; ++it )
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
    
    bool geometryTypeValid( const GeometryType &type ) const
    {
      int codim = dimension - type.dim();
      const std :: vector< GeometryType > &gT = geomTypes( codim );
      for( size_t i = 0; i < gT.size(); ++i ) {
        if( gT[ i ] == type )
          return true;
      }
      return false;
    }
  };


  
  template< class GridType >
  class WrappedPeriodicLeafIndexSet
  : public IndexSetWrapper< PeriodicLeafIndexSet< GridType > >
  {
    // Marked for revision in grid/common/defaultindexsets.hh
    // Whatever it is good for!
    enum { myType = 5 };
    
    typedef PeriodicLeafIndexSet< GridType > IndexSetType;
      
  private:
    typedef WrappedPeriodicLeafIndexSet< GridType > ThisType;
    typedef IndexSetWrapper< IndexSetType > BaseType;

  protected:
    IndexSetType indexSet_;

  public:
    enum { ncodim = GridType :: dimension + 1 };

    //! Constructor for wrapper and wrapped object.
    WrappedPeriodicLeafIndexSet ( const GridType &grid )
    : BaseType( indexSet_ ),
      indexSet_( grid )
    {
    }

    inline const std :: string name () const
    {
      return indexSet_.name();
    }

    //! return type (for Grape In/Output)
    static int type ()
    {
      return myType;
    }

    /*! \brief returns reference to singleton
     *
     *  Do we get a problem, if called with two different grids?
     */
    static ThisType &instance ( const GridType &grid )
    {
      static ThisType set( grid );
      return set;
    }
  };

}

#endif
