#ifndef DUNE_FEM_SPACE_MAPPER_GHOST_HH
#define DUNE_FEM_SPACE_MAPPER_GHOST_HH

#include <cstddef>

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/space/mapper/capabilities.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class GridPart, class Mapper >
    class AuxiliaryDofs;

    template< class AuxiliaryDofs >
    struct PrimaryDofs;



    namespace __GhostDofMapper
    {

      // BuildDataHandle
      // ---------------

      template< class BaseMapper >
      struct BuildDataHandle
        : public CommDataHandleIF< BuildDataHandle< BaseMapper >, std::pair< int, std::size_t > >
      {
        typedef std::pair< int, std::size_t > Data;

        explicit BuildDataHandle ( int rank, const BaseMapper &baseMapper, std::vector< std::tuple< int, std::size_t, std::size_t > > &masters )
          : rank_( rank ), baseMapper_( baseMapper ), masters_( masters )
        {}

        bool contains ( int dim, int codim ) const { return baseMapper_.contains( codim ); }
        bool fixedSize ( int dim, int codim ) const { return baseMapper_.fixedDataSize( codim ); }

        template< class Buffer, class Entity >
        void gather ( Buffer &buffer, const Entity &entity ) const
        {
          baseMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto index ) {
              std::get< 0 >( masters_[ index ] ) = rank_;
              buffer.write( Data( std::get< 0 >( masters_[ index ] ), std::get< 1 >( masters_[ index ] ) ) );
            } );
        }

        template< class Buffer, class Entity >
        void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
        {
          assert( n == size( entity ) );

          baseMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto index ) {
              Data remote( -1, index );
              buffer.read( remote );
              assert( remote.first >= 0 );

              auto &local = masters_[ index ];
              if( (std::get< 0 >( local ) < 0) || (remote.first < std::get< 0 >( local )) )
                std::tie( std::get< 0 >( local ), std::get< 1 >( local ) ) = remote;
          } );
        }

        template< class Entity >
        std::size_t size ( const Entity &entity ) const
        {
          return baseMapper_.numEntityDofs( entity );
        }

      protected:
        int rank_;
        const BaseMapper &baseMapper_;
        std::vector< std::tuple< int, std::size_t, std::size_t > > &masters_;
      };



      // ConstIterator
      // -------------

      template< class Index >
      class ConstIterator
      {
        typedef ConstIterator< Index > ThisType;

      public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef Index value_type;
        typedef Index difference_type;
        typedef Envelope< Index > pointer;
        typedef Index reference;

        ConstIterator () noexcept = default;
        explicit ConstIterator ( Index index ) noexcept : index_( index ) {}

        Index operator* () const noexcept { return index_; }
        Envelope< Index > operator-> () const noexcept { return Envelope< Index >( index_ ); }

        Index operator[] ( Index n ) const noexcept { return index_ + n; }

        bool operator== ( const ThisType &other ) const noexcept { return (index_ == other.index_); }
        bool operator!= ( const ThisType &other ) const noexcept { return (index_ != other.index_); }

        ThisType &operator++ () noexcept { ++index_; return *this; }
        ThisType operator++ ( int ) noexcept { ThisType copy( *this ); ++(*this); return copy; }

        ThisType &operator-- () noexcept { --index_; return *this; }
        ThisType operator-- ( int ) noexcept { ThisType copy( *this ); --(*this); return copy; }

        ThisType &operator+= ( Index n ) noexcept { index_ += n; return *this; }
        ThisType &operator-= ( Index n ) noexcept { index_ -= n; return *this; }

        ThisType operator+ ( Index n ) const noexcept { return ThisType( index_ + n ); }
        ThisType operator- ( Index n ) const noexcept { return ThisType( index_ - n ); }

        friend ThisType operator+ ( Index n, const ThisType &i ) noexcept { return i + n; }

        Index operator- ( const ThisType &other ) const noexcept { return (index_ - other.index_); }

        bool operator< ( const ThisType &other ) const noexcept { return (index_ < other.index_); }
        bool operator<= ( const ThisType &other ) const noexcept { return (index_ <= other.index_); }
        bool operator>= ( const ThisType &other ) const noexcept { return (index_ >= other.index_); }
        bool operator> ( const ThisType &other ) const noexcept { return (index_ > other.index_); }

      private:
        Index index_ = 0;
      };

    } // namespace __GhostDofMapper



    // GhostDofMapper
    // --------------

    template< class GridPart, class BaseMapper, class GlobalKey = std::size_t >
    class GhostDofMapper
    {
      typedef GhostDofMapper< GridPart, BaseMapper, GlobalKey > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef BaseMapper BaseMapperType;

      typedef std::size_t SizeType;
      typedef GlobalKey GlobalKeyType;

      typedef typename BaseMapperType::ElementType ElementType;

      GhostDofMapper ( const GridPartType &gridPart, BaseMapperType &baseMapper )
        : gridPart_( gridPart ), baseMapper_( baseMapper )
      {
        update();
      }

      GhostDofMapper ( const ThisType & ) = delete;
      GhostDofMapper ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        baseMapper().mapEach( element, [ this, &f ] ( auto local, auto i ) { f( local, mapping_[ i ] ); } );
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
        baseMapper().mapEachEntityDof( entity, [ this, &f ] ( auto local, auto i ) { f( local, mapping_[ i ] ); } );
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

      SizeType interiorSize () const { return interiorSize_; }
      SizeType ghostSize () const { return ghostSize_; }

      SizeType size () const { return interiorSize() + ghostSize(); }

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
        std::size_t baseSize = baseMapper().size();
        mapping_.resize( baseSize );

        std::vector< std::tuple< int, std::size_t, std::size_t > > masters( baseSize );
        for( std::size_t i = 0; i < baseSize; ++i )
          masters[ i ] = std::make_tuple( -1, i, i );

        const int rank = gridPart().comm().rank();
        __GhostDofMapper::BuildDataHandle< BaseMapper > dataHandle( rank, baseMapper_, masters );
        gridPart().communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
        // at this point all shared DoFs are assigned to their master rank
        // all other DoFs (with rank -1) are not shared at all

        // assign indices to interiors
        interiorSize_ = ghostSize_ = 0;
        for( const auto &m : masters )
        {
          if( (std::get< 0 >( m ) == -1) || (std::get< 0 >( m ) == rank ) )
            mapping_[ std::get< 2 >( m ) ] = interiorSize_++;
          else
            masters[ ghostSize_++ ] = m;
        }
        masters.resize( ghostSize_ );

        // sort the masters (by the first two components) to find duplicate ghosts
        const auto less = [] ( auto a, auto b ) { return (std::get< 0 >( a ) < std::get< 0 >( b )) || ((std::get< 0 >( a ) == std::get< 0 >( b )) && (std::get< 1 >( a ) < std::get< 1 >( b ))); };
        std::sort( masters.begin(), masters.end(), less );

        // assign indices to ghosts
        ghostSize_ = 0;
        std::tuple< int, std::size_t, std::size_t > current( -1, 0, 0 );
        for( const auto &m : masters )
        {
          if( less( current, m ) )
          {
            current = m;
            std::get< 2 >( current ) = interiorSize_ + ghostSize_++;
          }
          mapping_[ std::get< 2 >( m ) ] = std::get< 2 >( current );
        }
      }

      const GridPartType &gridPart () const { return gridPart_; }
      const BaseMapperType &baseMapper () const { return baseMapper_; }

      const std::vector< GlobalKeyType > &mapping () const { return mapping_; }

    private:
      const GridPartType &gridPart_;
      BaseMapperType &baseMapper_;
      std::vector< GlobalKeyType > mapping_;
      SizeType interiorSize_, ghostSize_;
    };



    // Capabilities for IndexSetDofMapper
    // ----------------------------------

    namespace Capabilities
    {

      template< class GridPart, class BaseMapper, class GlobalKey >
      struct isAdaptiveDofMapper< GhostDofMapper< GridPart, BaseMapper, GlobalKey > >
      {
        static const bool v = false;
      };

      template< class GridPart, class BaseMapper, class GlobalKey >
      struct isConsecutiveIndexSet< GhostDofMapper< GridPart, BaseMapper, GlobalKey > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities



    // AuxiliaryDofs for GhostDofMapper
    // ----------------------------

    template< class GridPart, class BaseMapper, class GlobalKey >
    class AuxiliaryDofs< GridPart, GhostDofMapper< GridPart, BaseMapper, GlobalKey > >
    {
      typedef AuxiliaryDofs< GridPart, GhostDofMapper< GridPart, BaseMapper, GlobalKey > > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef GhostDofMapper< GridPart, BaseMapper, GlobalKey > MapperType;

      typedef typename MapperType::GlobalKeyType GlobalKeyType;
      typedef typename MapperType::SizeType SizeType;

      typedef __GhostDofMapper::ConstIterator< GlobalKeyType > ConstIteratorType;

      explicit AuxiliaryDofs ( const MapperType &mapper )
        : mapper_( mapper )
      {}

      AuxiliaryDofs ( const GridPartType &gridPart, const MapperType &mapper )
        : AuxiliaryDofs( mapper )
      {}

      //! return dof number of salve with index
      GlobalKeyType operator [] ( int index ) const { return mapper().interiorSize() + index; }

      //! return number of auxiliary dofs
      SizeType size () const { return mapper().ghostSize()+1; }

      //! return number of primary dofs
      SizeType primarySize () const { return mapper().interiorSize(); }

      ConstIteratorType begin () const { return ConstIteratorType( mapper().interiorSize() ); }
      ConstIteratorType end () const { return ConstIteratorType( mapper().interiorSize() + mapper().ghostSize() ); }

      //! return true if index is contained, meaning it is an auxiliary dof
      bool contains ( GlobalKeyType index ) const { return (static_cast< SizeType >( index ) >= mapper().interiorSize()); }

      void rebuild () {}

      const MapperType &mapper () const { return mapper_; }
      const GridPartType &gridPart () const { return mapper().gridPart(); }

    private:
      const MapperType &mapper_;
    };


    /** \brief @ingroup Communication
     *
     *  Apply action encoded in Functor f to all primary dofs.
     *  \param[in]  auxiliaryDofs  AuxiliaryDofs instance (from space)
     *  \param[in]  f  a Functor or lambda offering operator()( const size_t dof )
     *
     *  \note: This is a specialization for the Petsc  GhostMapper based AuxiliaryDofs.
     */
    template <class GridPart, class BaseMapper, class GlobalKey, class F>
    static void forEachPrimaryDof( const AuxiliaryDofs< GridPart, GhostDofMapper< GridPart, BaseMapper, GlobalKey > >& auxiliaryDofs, F&& f )
    {
      const size_t size = auxiliaryDofs.mapper().interiorSize();
      for( size_t dof = 0 ; dof<size; ++dof )
      {
        // apply action to primary dof
        f( dof );
      }
    }




    // PrimaryDofs for AuxiliaryDofs< GhostDofMapper >
    // ------------------------------------------

    template< class GridPart, class BaseMapper, class GlobalKey >
    struct PrimaryDofs< AuxiliaryDofs< GridPart, GhostDofMapper< GridPart, BaseMapper, GlobalKey > > >
    {
      typedef AuxiliaryDofs< GridPart, GhostDofMapper< GridPart, BaseMapper, GlobalKey > > AuxiliaryDofsType;

      typedef typename AuxiliaryDofsType::GlobalKeyType GlobalKeyType;
      typedef typename AuxiliaryDofsType::GridPartType GridPartType;
      typedef typename AuxiliaryDofsType::MapperType MapperType;
      typedef typename AuxiliaryDofsType::SizeType SizeType;

      typedef __GhostDofMapper::ConstIterator< GlobalKeyType > ConstIteratorType;

      [[deprecated("Use forEachPrimaryDof instead!")]]
      explicit PrimaryDofs ( const AuxiliaryDofsType &auxiliaryDofs )
        : mapper_( auxiliaryDofs.mapper() )
      {}

      ConstIteratorType begin () const { return ConstIteratorType( 0 ); }
      ConstIteratorType end () const { return ConstIteratorType( size() ); }

      SizeType size () const { return mapper().interiorSize(); }

      const MapperType &mapper () const { return mapper_; }
      const GridPartType &gridPart () const { return mapper().gridPart(); }

    private:
      const MapperType &mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_MAPPER_GHOST_HH
