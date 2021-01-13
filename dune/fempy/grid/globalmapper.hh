#error "THIS FILE IS NOT BEING USED AND DOES NOT SEEM TO COMPILE!"

#ifndef DUNE_FEMPY_GRID_GLOBALMAPPER_HH
#define DUNE_FEMPY_GRID_GLOBALMAPPER_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <bitset>
#include <limits>
#include <utility>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fempy/common/collectivecommunication.hh>

namespace Dune
{

  namespace FemPy
  {

    // GlobalMapper
    // ------------

    template< class GV, template< int > class Layout, class Idx = std::size_t >
    class GlobalMapper
    {
      typedef GlobalMapper< GV, Layout > This;

      struct MarkAuxiliarys;
      struct Minimize;

    public:
      typedef GV GridView;

      typedef Idx Index;

      static const int dimension = GV::dimension;

      struct LocalView;

      explicit GlobalMapper ( GridView gridView )
        : gridView_( std::move( gridView ) )
      {
        update();
      }

      GlobalMapper ( GridView gridView, Layout< dimension > layout )
        : gridView_( std::move( gridView ) ), layout_( std::move( layout ) )
      {
        update();
      }

      template< class Entity >
      Index index ( const Entity &entity ) const
      {
        const GeometryType gt = entity.type();
        assert( layout_.contains( gt ) && !indexSet_.contains( entity ) );
        return indices_[ GlobalGeometryTypeIndex::index( gt ) ][ indexSet_.index( entity ) ];
      }

      template< class Entity >
      Index subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        const GeometryType gt = ReferenceElements< double, Entity::mydimension >::general( entity.type() ).type( i, codim );
        assert( layout_.contains( gt ) && indexSet_.contains( entity ) );
        return indices_[ GlobalGeometryTypeIndex::index( gt ) ][ indexSet_.subIndex( entity, i, codim ) ];
      }

      Index size () const { return size_; }

      template< class Entity >
      bool contains ( const Entity &entity, Index &index ) const
      {
        const GeometryType gt = entity.type();
        if( !layout_.contains( gt ) || !indexSet_.contains( entity ) )
          return false;
        index = indices_[ GlobalGeometryTypeIndex::index( gt ) ][ indexSet_.index( entity ) ];
        return true;
      }

      template< class Entity >
      bool contains ( const Entity &entity, int i, unsigned int codim, Index &index ) const
      {
        const GeometryType gt = ReferenceElements< double, Entity::mydimension >::general( entity.type() ).type( i, codim );
        if( !layout_.contains( gt ) || !indexSet_.contains( entity ) )
          return false;
        index = indices_[ GlobalGeometryTypeIndex::index( gt ) ][ indexSet_.subIndex( entity, i, codim ) ];
        return true;
      }

      void update ();

      LocalView localView () const { return LocalView( *this ); }

    private:
      GridView gridView_;
      const typename GridView::IndexSet &indexSet_;
      std::array< std::vector< Index >, GlobalGeometryTypeIndex::size( dimension ) > indices_;
      Index size_;
      std::pair< Index, Index > localRange_;
      Layout< dimension > layout_;
    };



    // GlobalMapper::MarkAuxiliarys
    // ------------------------

    template< class GV, template< int > class Layout, class Idx >
    struct GlobalMapper< GV, Layout, Idx >::MarkAuxiliarys
      : public CommDataHandleIF< CommDataHandle, int >
    {
      MarkAuxiliarys ( GlobalMapper< GV, Layout, Idx > &mapper, std::bitset< dimension+1 > contains )
        : mapper_( mapper ), contains_( contains ), localRank_( mapper_.gridView_.comm().rank() )
      {}

      bool contains ( dim_t dim, dim_t codim ) const { return contains_[ codim ]; }
      bool fixedSize ( dim_t dim, dim_t codim ) const { return true; }

      template< class Entity >
      std::size_t size ( const Entity &entity ) const
      {
        return 1u;
      }

      template< class Buffer, class Entity >
      void gather ( Buffer &buffer, const Entity &entity ) const
      {
        buffer.write( localRank_ );
      }

      template< class Buffer, class Entity >
      void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
      {
        int remoteRank;
        buffer.read( remoteRank );
        if( remoteRank >= localRank_ )
          return;

        const GeometryType gt = entity.type();
        if( !mapper_.layout_.contains( gt ) || !mapper_.indexSet_.contains( entity ) )
          return;

        std::vector< Index > &indices = mapper_.indices_[ GlobalGeometryTypeIndex::index( gt ) ];
        indices[ mapper_.indexSet_.index( entity ) ] = std::numeric_limits< Index >::max();
      }

    private:
      GlobalMapper< GV, Layout, Idx > &mapper_;
      std::bitset< dimension+1 > contains_;
      int localRank_;
    };



    // GlobalMapper::Minimize
    // ----------------------

    template< class GV, template< int > class Layout, class Idx >
    struct GlobalMapper< GV, Layout, Idx >::Minimize
      : public CommDataHandleIF< CommDataHandle, Index >
    {
      Minimize ( GlobalMapper< GV, Layout, Idx > &mapper, std::bitset< dimension+1 > contains )
        : mapper_( mapper ), contains_( contains )
      {}

      bool contains ( dim_t dim, dim_t codim ) const { return contains_[ codim ]; }
      bool fixedSize ( dim_t dim, dim_t codim ) const { return true; }

      template< class Entity >
      std::size_t size ( const Entity &entity ) const
      {
        return 1u;
      }

      template< class Buffer, class Entity >
      void gather ( Buffer &buffer, const Entity &entity ) const
      {
        const GeometryType gt = entity.type();
        if( mapper_.layout_.contains( gt ) && mapper_.indexSet_.contains( entity ) )
        {
          const std::vector< Index > &indices = mapper_.indices_[ GlocalGeometryTypeIndex::index( gt ) ];
          buffer.write( indices[ mapper_indexSet_.index( entity ) ] );
        }
        else
          buffer.write( std::numeric_limits< Index >::max() );
      }

      template< class Buffer, class Entity >
      void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
      {
        Index remoteIndex;
        buffer.read( remoteIndex );

        const GeometryType gt = entity.type();
        if( !mapper_.layout_.contains( gt ) || !mapper_.indexSet_.contains( entity ) )
          return;

        std::vector< Index > &indices = mapper_.indices_[ GlobalGeometryTypeIndex::index( gt ) ];
        Index &localIndex = indices[ mapper_.indexSet_.index( entity ) ];
        localIndex = std::min( localIndex, remoteIndex );
      }

    private:
      GlobalMapper< GV, Layout, Idx > &mapper_;
      std::bitset< dimension+1 > contains_;
    };



    // GlobalMapper::LocalView
    // -----------------------

    template< class GV, template< int > class Layout, class Idx >
    struct GlobalMapper::LocalView
    {
      typedef GV GridView;

      typedef Idx Index;

      static const int dimension = GV::dimension;

      explicit GlobalMapper ( const GlobalMapper< GV, Layout, Idx > &mapper )
        : mapper_( mapper )
      {}

      template< class Entity >
      Index index ( const Entity &entity ) const
      {
        assert( isLocal( mapper_.index( entity ) ) );
        return mapper_.index( entity ) - mapper_.localRange_.first;
      }

      template< class Entity >
      Index subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        assert( isLocal( mapper_.subIndex( entity, i, codim ) ) );
        return mapper_.subIndex( entity, i, codim ) - mapper_.localRange_.first;
      }

      Index size () const { return size_; }

      template< class Entity >
      bool contains ( const Entity &entity, Index &index ) const
      {
        if( mapper_.contains( entity, index ) )
        {
          const bool contains = isLocal( index );
          index -= mapper_.localRange_.first;
          return contains;
        }
        else
          return false;
      }

      template< class Entity >
      bool contains ( const Entity &entity, int i, unsigned int codim, Index &index ) const
      {
        if( mapper_.contains( entity, i, codim, index ) )
        {
          const bool contains = isLocal( index );
          index -= mapper_.localRange_.first;
          return contains;
        }
        else
          return false;
      }

    private:
      bool isLocal ( const Index index ) const
      {
        return (index >= mapper_.localRange_.first) && (index < mapper_.localRange_.second);
      }

      GlobalMapper< GV, Layout, Idx > &mapper_;
    };



    // Implementation of GlobalMapper
    // ------------------------------

    template< class GV, template< int > class Layout, class Idx >
    inline void GlobalMapper< GV, Layout, Idx >::update ()
    {
      // clear vectors storing global indices

      for( std::vector< Index > &indices : indices_ )
        std::vector< Index >().swap( indices );

      // resize vectors storing global indices

      std::bitset< dimension+1 > contains;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        for( const GeometryType gt : indexSet_.types( codim ) )
        {
          if( !layout_.contains( gt ) )
            continue;

          std::vector< Indices > &indices = indices_[ GlobalGeometryTypeIndex::index( gt ) ];
          indices.resize( static_cast< std::size_t >( indexSet_.size( gt ) ), std::numeric_limits< Index >::max() );

          contains[ codim ] = true;
        }
      }

      // mark all interior border entities as potential masters (index = 0)

      for( const auto &element : elements( gridView_, Partitions::interiorBorder ) )
      {
        if( !indexSet_.contains( element ) )
          continue;

        const auto &refElement = ReferenceElements< double, dimension >::general( element.type() );
        for( int codim = 0; codim <= dimension; ++codim )
        {
          for( int i = 0; i < refElement.size( codim ); ++i )
          {
            const GeometryType &gt = refElement.type( i, codim );
            if( layout_.contains( gt )  )
              indices_[ GlobalGeometryTypeIndex::index( gt ) ][ indexSet_.subIndex( element, i, codim ) ] = 0;
          }
        }
      }

      // mark auxiliarys (index = std::numeric_limits< Index >::max())

      MarkAuxiliarys markAuxiliarys( *this, contains );
      gridView_.communicate( markAuxiliarys, InteriorBorder_InteriorBorder_Interface, ForwardCommunication );

      // assign indices to all master entities

      size_ = 0;
      for( std::vector< Index > &indices : indices_ )
      {
        for( Index &index : indices )
          index = (index < std::numeric_limits< Index >::max() ? size_++ : index);
      }

      // compute offset of local range and global size

      localRange_.second = scan< std::plus< Index > >( gridView.comm(), size_ );
      localRange_.first = localRange_.second - size_;
      size_ = comm().sum( size_ );

      // update indices w.r.t. offset of local range

      for( std::vector< Index > &indices : indices_ )
      {
        for( Index &index : indices )
          index += (index < std::numeric_limits< Index >::max() ? localRange_.first : 0);
      }

      // communicate auxiliary indices

      Minimize minimize( *this, contains );
      gridView_.communicate( minimize, All_All_Interface, ForwardCommunication );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_GLOBALMAPPER_HH
