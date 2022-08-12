#ifndef DUNE_FEM_SPACE_MAPPER_PARALLEL_HH
#define DUNE_FEM_SPACE_MAPPER_PARALLEL_HH

#include <cstddef>

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/space/common/auxiliarydofs.hh>
#include <dune/fem/space/mapper/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    namespace __ParallelDofMapper
    {

      // BuildDataHandle
      // ---------------

      template< class GridPart, class BaseMapper, class GlobalKey >
      struct BuildDataHandle
        : public CommDataHandleIF< BuildDataHandle< GridPart, BaseMapper, GlobalKey >, GlobalKey >
      {
        explicit BuildDataHandle ( const BaseMapper &baseMapper, const AuxiliaryDofs< GridPart, BaseMapper > &auxiliaryDofs, std::vector< GlobalKey > &mapping )
          : baseMapper_( baseMapper ), auxiliaryDofs_( auxiliaryDofs ), mapping_( mapping )
        {}

        bool contains ( int dim, int codim ) const { return baseMapper_.contains( codim ); }
        bool fixedSize ( int dim, int codim ) const { return false; }

        template< class Buffer, class Entity >
        void gather ( Buffer &buffer, const Entity &entity ) const
        {
          baseMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto index ) {
              if( !auxiliaryDofs_.contains( index ) )
                buffer.write( mapping_[ index ] );
            } );
        }

        template< class Buffer, class Entity >
        void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
        {
          if( n == 0 )
            return;

          assert( n == static_cast< std::size_t >( baseMapper_.numEntityDofs( entity ) ) );
          baseMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto index ) {
              assert( auxiliaryDofs_.contains( index ) );
              buffer.read( mapping_[ index ] );
          } );
        }

        template< class Entity >
        std::size_t size ( const Entity &entity ) const
        {
          std::size_t size = 0;
          baseMapper_.mapEachEntityDof( entity, [ this, &size ] ( int, auto index )
              { size += static_cast< std::size_t >( !auxiliaryDofs_.contains( index ) ); } );
          return size;
        }

      protected:
        const BaseMapper &baseMapper_;
        const AuxiliaryDofs< GridPart, BaseMapper > &auxiliaryDofs_;
        std::vector< GlobalKey > &mapping_;
      };

    } // namespace __ParallelDofMapper



    // ParallelDofMapper
    // -----------------

    template< class GridPart, class BaseMapper, class GlobalKey = std::size_t >
    class ParallelDofMapper
    {
      typedef ParallelDofMapper< GridPart, BaseMapper, GlobalKey > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef BaseMapper BaseMapperType;

      typedef std::size_t SizeType;
      typedef GlobalKey GlobalKeyType;

      typedef typename BaseMapperType::ElementType ElementType;

      ParallelDofMapper ( const GridPartType &gridPart, const BaseMapperType &baseMapper )
        : gridPart_( gridPart ), baseMapper_( baseMapper )
      {
        update();
      }

      ParallelDofMapper ( const ThisType & ) = delete;
      ParallelDofMapper ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        baseMapper().mapEach( element, [ this, f ] ( auto local, auto i ) { f( local, mapping_[ i ] ); } );
      }

      void map ( const ElementType &element, std::vector< GlobalKeyType > &indices ) const
      {
        indices.resize( numDofs( element ) );
        mapEach( element, [ &indices ] ( int local, GlobalKeyType global ) { indices[ local ] = global; } );
      }

      void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
      {
        baseMapper().onSubEntity( element, i, c, indices );
      }

      unsigned int maxNumDofs () const { return baseMapper().maxNumDofs(); }
      unsigned int numDofs ( const ElementType &element ) const { return baseMapper().numDofs( element ); }

      // assignment of DoFs to entities

      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        baseMapper().mapEachEntityDof( entity, [ this, f ] ( auto local, auto i ) { f( local, mapping_[ i ] ); } );
      }

      template< class Entity >
      void mapEntityDofs ( const Entity &entity, std::vector< GlobalKeyType > &indices ) const
      {
        indices.resize( numEntityDofs( entity ) );
        mapEachEntityDof( entity, [ &indices ] ( int local, GlobalKeyType global ) { indices[ local ] = global; } );
      }

      template< class Entity >
      unsigned int numEntityDofs ( const Entity &entity ) const
      {
        return baseMapper().numEntityDofs( entity );
      }

      // global information

      bool contains ( int codim ) const { return baseMapper().contains( codim ); }

      bool fixedDataSize ( int codim ) const { return baseMapper().fixedDataSize( codim ); }

      SizeType size () const { return size_; }

      // adaptation interface

      bool consecutive () const { return false; }

      int numBlocks () const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }
      SizeType offSet ( int blk ) const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }
      SizeType oldOffSet ( int blk ) const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }
      SizeType numberOfHoles ( int blk ) const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }
      SizeType oldIndex ( SizeType hole, int blk ) const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }
      SizeType newIndex ( SizeType hole, int blk ) const { DUNE_THROW( NotImplemented, "Adaptive dof mapper interface not implemented." ); }

      // update

      void update ()
      {
        AuxiliaryDofs< GridPartType, BaseMapperType > auxiliaryDofs( gridPart(), baseMapper() );
        auxiliaryDofs.rebuild();
        //auto primaryDofs = Dune::Fem::primaryDofs( auxiliaryDofs );

        const auto primarySize = auxiliaryDofs.primarySize(); // primaryDofs.size();
        size_ = primarySize;  // primaryDofs.size();
        offset_ = exScan( gridPart().comm(), size_ );
        size_ = gridPart().comm().sum( size_ );

        std::size_t baseSize = baseMapper().size();
        mapping_.resize( baseSize );
        GlobalKeyType next = static_cast< GlobalKeyType >( offset_ );
        std::vector< GlobalKeyType >& mapping = mapping_;
        auto mapNext = [&mapping, &next] (const auto i) { mapping[ i ] = next++; };
        // for all primary dofs build mapping
        forEachPrimaryDof( auxiliaryDofs, mapNext );
        //for( const auto i : primaryDofs )
        //  mapping_[ i ] = next++;
        assert( next == static_cast< GlobalKeyType >( offset_ + primarySize ) );

        __ParallelDofMapper::BuildDataHandle< GridPartType, BaseMapperType, GlobalKeyType > dataHandle( baseMapper(), auxiliaryDofs, mapping_ );
        gridPart().communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
      }

      const GridPartType &gridPart () const { return gridPart_; }
      const BaseMapperType &baseMapper () const { return baseMapper_; }

      const std::vector< GlobalKeyType > &mapping () const { return mapping_; }

    private:
      template< class Comm, class T >
      static T exScan ( const Communication< Comm > &comm, T in )
      {
        return T( 0 );
      }

#if HAVE_MPI
      template< class T >
      static T exScan ( const Communication< MPI_Comm > &comm, T in )
      {
        T out( 0 );
        MPI_Exscan( &in, &out, 1, MPITraits< T >::getType(), MPI_SUM, static_cast< MPI_Comm >( comm ) );
        return out;
      }
#endif // #if HAVE_MPI

      const GridPartType &gridPart_;
      const BaseMapperType &baseMapper_;
      std::vector< GlobalKeyType > mapping_;
      SizeType offset_, size_;
    };



    // Capabilities for IndexSetDofMapper
    // ----------------------------------

    namespace Capabilities
    {

      template< class GridPart, class BaseMapper, class GlobalKey >
      struct isAdaptiveDofMapper< ParallelDofMapper< GridPart, BaseMapper, GlobalKey > >
      {
        static const bool v = false;
      };

      template< class GridPart, class BaseMapper, class GlobalKey >
      struct isConsecutiveIndexSet< ParallelDofMapper< GridPart, BaseMapper, GlobalKey > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_MAPPER_PARALLEL_HH
