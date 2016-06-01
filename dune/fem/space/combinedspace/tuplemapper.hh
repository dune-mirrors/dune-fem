#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH

#include <array>
#include <tuple>
#include <utility>

#include <dune/common/forloop.hh>
#include <dune/common/tuples.hh>

#include <dune/common/std/utility.hh>
#include <dune/fem/common/utility.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>


namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------
    template< class GridPart, class ... Mapper >
    class TupleMapper;


#ifndef DOXYGEN

    namespace __TupleMapper
    {

      // Traits
      // ------

      template< class GridPart, class ... Mapper >
      struct Traits
      {
        static_assert( Std::are_all_same< typename Mapper::ElementType ... >::value,
                       "TupleMapper needs common ElementType" );

        typedef typename std::tuple_element< 0, std::tuple< Mapper ... > >::type FirstMapperType;
        typedef typename FirstMapperType::ElementType ElementType;
        typedef typename FirstMapperType::SizeType SizeType;
        typedef typename FirstMapperType::GlobalKeyType GlobalKeyType;

        typedef TupleMapper< GridPart, Mapper ... > DofMapperType;
      };

      // CombinedIndex
      // -------------

      template< class Index, class Int, Int i >
      struct CombinedIndex
      {
        constexpr CombinedIndex ( Index index, Index offset ) : index_( index ), offset_( offset ) {}

        static constexpr Int component () { return i; }

        constexpr operator Index () const { return index_ + offset_; }

        constexpr Index index () const { return index_; }

        constexpr Index offset () const { return offset_; }

      private:
        Index index_, offset_;
      };


      // DofMapper
      // ---------

      template< class T, template< class > class Base = Dune::Fem::DofMapper >
      class DofMapper;

      template< class GridPart, class ... Mapper, template< class > class Base >
      class DofMapper< Traits< GridPart, Mapper ... >, Base >
        : public Base< Traits< GridPart, Mapper ... > >
      {
        typedef Base< Traits< GridPart, Mapper ... > > BaseType;

        // helper classes
        template< int >
        struct ComputeOffSet;
        template< int >
        struct MapEach;
        template< int >
        struct MapEachEntityDof;

        // FunctorWrapper
        // --------------

        template< class Functor, int i >
        struct FunctorWrapper
        {
          FunctorWrapper ( Functor functor, int localOffset, int globalOffset )
            : functor_( functor ),
              localOffset_( localOffset ),
              globalOffset_( globalOffset )
          {}

          template< class GlobalKey >
          void operator() ( int localDof, const GlobalKey &globalKey )
          {
            functor_( localDof + localOffset_, CombinedIndex< GlobalKey, int, i >( globalKey, globalOffset_ ) );
          }

          template< class GlobalKey >
          void operator() ( const GlobalKey &globalKey )
          {
            functor_( CombinedIndex< GlobalKey, int, i >( globalKey, globalOffset_ ) );
          }

        private:
          Functor functor_;
          const int localOffset_;
          const int globalOffset_;
        };

        // size of the Mapper Tuple
        static const int mapperTupleSize = sizeof ... ( Mapper );

        typedef std::array< typename BaseType::SizeType, mapperTupleSize + 1 > OffsetType;

      public:
        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;
        typedef typename BaseType::Traits::GlobalKeyType GlobalKeyType;

        typedef GridPart GridPartType;

        DofMapper ( GridPartType &gridPart, Mapper & ... mapper )
          : gridPart_( gridPart ),
            mapperTuple_( mapper ... )
        {
          init();
        }

        DofMapper ( GridPartType &gridPart, Mapper && ... mapper )
          : gridPart_( gridPart ),
            mapperTuple_( std::move( mapper ) ... )
        {
          init();
        }

        SizeType size () const { return size( Std::index_sequence_for< Mapper ... >() ); }

        bool contains ( const int codim ) const { return contains( codim, Std::index_sequence_for< Mapper ... >() ); }

        bool fixedDataSize ( int codim ) const { return fixedDataSize( codim, Std::index_sequence_for< Mapper ... >() ); }

        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const
        {
          OffsetType localOffset;
          localOffset[ 0 ] = 0;
          ForLoop< MapEach, 0, mapperTupleSize - 1 >::apply( localOffset, globalOffset_, element, f, mapperTuple_ );
        }

        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const
        {
          OffsetType localOffset;
          localOffset[ 0 ] = 0;
          ForLoop< MapEachEntityDof, 0, mapperTupleSize - 1 >::apply( localOffset, globalOffset_, entity, f, mapperTuple_ );
        }

        int maxNumDofs () const { return maxNumDofs( Std::index_sequence_for< Mapper ... >() ); }

        SizeType numDofs ( const ElementType &element ) const { return numDofs( element, Std::index_sequence_for< Mapper ... >() ); }

        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const { return numEntityDofs( entity, Std::index_sequence_for< Mapper ... >() ); }


        static constexpr bool consecutive () noexcept { return false; }

        SizeType numBlocks () const
        {
          DUNE_THROW( NotImplemented, "Method numBlocks() called on non-adaptive block mapper" );
        }

        SizeType numberOfHoles ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method numberOfHoles() called on non-adaptive block mapper" );
        }

        GlobalKeyType oldIndex ( int hole, int ) const
        {
          DUNE_THROW( NotImplemented, "Method oldIndex() called on non-adaptive block mapper" );
        }

        GlobalKeyType newIndex ( int hole, int ) const
        {
          DUNE_THROW( NotImplemented, "Method newIndex() called on non-adaptive block mapper" );
        }

        SizeType oldOffSet ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method oldOffSet() called on non-adaptive block mapper" );
        }

        SizeType offSet ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method offSet() called on non-adaptive block mapper" );
        }

        /*** NonInterface Methods ***/

        SizeType offset ( int i ) const { return globalOffset_[ i ]; }

        template< int i >
        SizeType subSize () const { return std::get< i >( mapperTuple_ ).size(); }

      protected:
        template< std::size_t ... i >
        SizeType size ( Std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).size() ... );
        }

        template< std::size_t ... i >
        bool fixedDataSize ( const int codim, Std::index_sequence< i ... > ) const
        {
          return Std::And( std::get< i >( mapperTuple_ ).fixedDataSize( codim ) ... );
        }

        template< std::size_t ... i >
        bool contains ( const int codim, Std::index_sequence< i ... > ) const
        {
          return Std::Or( std::get< i >( mapperTuple_ ).contains( codim ) ... );
        }

        template< std::size_t ... i >
        int maxNumDofs ( Std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).maxNumDofs() ... );
        }

        template< std::size_t ... i >
        SizeType numDofs ( const ElementType &element, Std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numDofs( element ) ... );
        }

        template< class Entity, std::size_t ... i >
        SizeType numEntityDofs ( const Entity &entity, Std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numEntityDofs( entity ) ... );
        }

        void init ()
        {
          globalOffset_[ 0 ] = 0;
          // compute new offsets
          ForLoop< ComputeOffSet, 0, mapperTupleSize - 1 >::apply( globalOffset_, mapperTuple_ );
        }

        GridPartType &gridPart_;
        std::tuple< Mapper ... > mapperTuple_;
        OffsetType globalOffset_;
      };


      // ComputeOffSet
      // -------------

      template< class GridPart, class ... Mapper, template< class > class Base >
      template< int i >
      struct DofMapper< Traits< GridPart, Mapper ... >, Base >::
      ComputeOffSet
      {
        template< class Tuple >
        static void apply ( OffsetType &offset, const Tuple &tuple )
        {
          offset[ i + 1 ] = offset[ i ] + std::get< i >( tuple ).size();
        }
      };


      // MapEachEntityDof
      // ----------------

      template< class GridPart, class ... Mapper, template< class > class Base >
      template< int i >
      struct DofMapper< Traits< GridPart, Mapper ... >, Base >::
      MapEachEntityDof
      {
        template< class Entity, class Functor, class Tuple >
        static void apply ( OffsetType &localOffset, const OffsetType &globalOffset, const Entity &entity, Functor f, const Tuple &tuple )
        {
          FunctorWrapper< Functor, i > wrappedFunctor( f, localOffset[ i ], globalOffset[ i ] );
          std::get< i >( tuple ).mapEachEntityDof( entity, wrappedFunctor );
          localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( tuple ).numEntityDofs( entity );
        }
      };


      // MapEach
      // -------

      template< class GridPart, class ... Mapper, template< class > class Base >
      template< int i >
      struct DofMapper< Traits< GridPart, Mapper ... >, Base >::
      MapEach
      {
        template< class Functor, class Tuple >
        static void apply ( OffsetType &localOffset, const OffsetType &globalOffset, const ElementType &element, Functor f, const Tuple &tuple )
        {
          FunctorWrapper< Functor, i > wrappedFunctor( f, localOffset[ i ], globalOffset[ i ] );
          std::get< i >( tuple ).mapEach( element, wrappedFunctor );
          localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( tuple ).numDofs( element );
        }
      };


      // AdaptiveDofMapper
      // -----------------

      template< class T >
      class AdaptiveDofMapper;

      template< class GridPart, class ... Mapper >
      class AdaptiveDofMapper< Traits< GridPart, Mapper ... > >
        : public DofMapper< Traits< GridPart, Mapper ... >, Dune::Fem::AdaptiveDofMapper >
      {
        typedef DofMapper< Traits< GridPart, Mapper ... >, Dune::Fem::AdaptiveDofMapper > BaseType;

      protected:

        // helper classes
        template< int >
        struct Offset;
        template< int >
        struct OldOffset;
        template< int >
        struct NewIndex;
        template< int >
        struct OldIndex;
        template< int >
        struct NumberOfHoles;

        // size of the Mapper Tuple
        static const int mapperTupleSize = sizeof ... ( Mapper );

        typedef std::array< typename BaseType::Traits::SizeType, mapperTupleSize + 1 > OffsetType;

        using BaseType::mapperTuple_;
        using BaseType::gridPart_;
        using BaseType::globalOffset_;

      public:
        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;
        typedef typename BaseType::GlobalKeyType GlobalKeyType;
        typedef GridPart GridPartType;

        AdaptiveDofMapper ( GridPartType &gridPart, Mapper & ... mapper )
          : BaseType( gridPart, mapper ... )
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        AdaptiveDofMapper ( GridPartType &gridPart, Mapper && ... mapper )
          : BaseType( gridPart, std::move( mapper ) ... )
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        ~AdaptiveDofMapper () { DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).removeIndexSet( *this ); }

        AdaptiveDofMapper ( const AdaptiveDofMapper & ) = delete;
        AdaptiveDofMapper ( AdaptiveDofMapper && ) = delete;

        static constexpr bool consecutive () noexcept { return true; }

        SizeType numBlocks () const { return numBlocks( Std::index_sequence_for< Mapper ... >() ); }

        SizeType numberOfHoles ( int block ) const
        {
          SizeType nHoles = 0;
          int comp = -1;
          ForLoop< NumberOfHoles, 0, mapperTupleSize - 1 >::apply( nHoles, comp, block, mapperTuple_ );
          return nHoles;
        }

        GlobalKeyType oldIndex ( int hole, int block ) const
        {
          int comp = -1;
          SizeType oIndex = 0;
          ForLoop< OldIndex, 0, mapperTupleSize - 1 >::apply( oIndex, comp, hole, block, mapperTuple_ );
          assert( comp >= 0 );
          return oIndex + globalOffset_[ comp ];
        }

        GlobalKeyType newIndex ( int hole, int block ) const
        {
          int comp = -1;
          SizeType nIndex = 0;
          ForLoop< NewIndex, 0, mapperTupleSize - 1 >::apply( nIndex, comp, hole, block, mapperTuple_ );
          assert( comp > 0 );
          return nIndex + globalOffset_[ comp ];
        }

        SizeType oldOffSet ( int block ) const
        {
          int comp = -1;
          SizeType oOffset = 0;
          ForLoop< OldOffset, 0, mapperTupleSize - 1 >::apply( oOffset, comp, block, mapperTuple_ );
          assert( comp >= 0 );
          return oOffset + oldGlobalOffset_[ comp ];
        }

        SizeType offSet ( int block ) const
        {
          int comp = -1;
          SizeType offset = 0;
          ForLoop< Offset, 0, mapperTupleSize - 1 >::apply( offset, comp, block, mapperTuple_ );
          assert( comp >= 0 );
          return offset + globalOffset_[ comp ];
        }


        void resize () { update(); }

        bool compress ()
        {
          update();
          return true;
        }

        void backup () const {}

        void restore () { update(); }

        template< class IOStream >
        void read ( IOStream &in ) { update(); }

        template< class IOStream >
        void write ( IOStream &out ) const {}

        template< class Entity >
        void insertEntity ( const Entity & ) { update(); }

        template< class Entity >
        void removeEntity ( const Entity & ) { update(); }

      protected:

        void update ()
        {
          oldGlobalOffset_ = globalOffset_;
          BaseType::init();
        }

        template< std::size_t ... i >
        SizeType numBlocks ( Std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numBlocks() ... );
        }

      private:
        OffsetType oldGlobalOffset_;
      };



      // NumberOfHoles
      // -------------

      template< class GridPart, class ... Mapper >
      template< int i >
      struct AdaptiveDofMapper< Traits< GridPart, Mapper ... > >::
      NumberOfHoles
      {
        template< class Tuple >
        static void apply ( SizeType &nHoles, int &comp, const int block, const Tuple &tuple )
        {
          if( comp >= 0 )
            return;
          const int localBlock = block - std::get< i >( tuple ).numBlocks();
          if( localBlock >= 0 )
          {
            comp = i;
            nHoles = std::get< i >( tuple ).numberOfHoles( localBlock );
          }
        }
      };



      // OldIndex
      // --------

      template< class GridPart, class ... Mapper >
      template< int i >
      struct AdaptiveDofMapper< Traits< GridPart, Mapper ... > >::
      OldIndex
      {
        template< class Tuple >
        static void apply ( SizeType &oldIndex, int &comp, const int hole, const int block, const Tuple &tuple )
        {
          if( comp >= 0 )
            return;
          const int localBlock = block - std::get< i >( tuple ).numBlocks();
          if( localBlock >= 0 )
          {
            comp = i;
            oldIndex = std::get< i >( tuple ).oldIndex( hole, localBlock );
          }
        }
      };


      // NewIndex
      // --------

      template< class GridPart, class ... Mapper >
      template< int i >
      struct AdaptiveDofMapper< Traits< GridPart, Mapper ... > >::
      NewIndex
      {
        template< class Tuple >
        static void apply ( SizeType &nIndex, int &comp, const int hole, const int block, const Tuple &tuple )
        {
          if( comp >= 0 )
            return;
          const int localBlock = block - std::get< i >( tuple ).numBlocks();
          if( localBlock >= 0 )
          {
            comp = i;
            nIndex = std::get< i >( tuple ).newIndex( hole, localBlock );
          }
        }
      };


      // OldOffset
      // ---------

      template< class GridPart, class ... Mapper >
      template< int i >
      struct AdaptiveDofMapper< Traits< GridPart, Mapper ... > >::
      OldOffset
      {
        template< class Tuple >
        static void apply ( SizeType &offset, int &comp, const int block, const Tuple &tuple )
        {
          if( comp >= 0 )
            return;
          const int localBlock = block - std::get< i >( tuple ).numBlocks();
          if( localBlock >= 0 )
          {
            comp = i;
            offset = std::get< i >( tuple ).oldOffSet( localBlock );
          }
        }
      };


      // Offset
      // ------

      template< class GridPart, class ... Mapper >
      template< int i >
      struct AdaptiveDofMapper< Traits< GridPart, Mapper ... > >::
      Offset
      {
        template< class Tuple >
        static void apply ( SizeType &offset, int &comp, const int block, const Tuple &tuple )
        {
          if( comp >= 0 )
            return;
          const int localBlock = block - std::get< i >( tuple ).numBlocks();
          if( localBlock >= 0 )
          {
            comp = i;
            offset = std::get< i >( tuple ).offSet( localBlock );
          }
        }
      };


      // Implementation
      // --------------

      template< class GridPart, class ... Mapper >
      struct Implementation
      {
        typedef typename std::conditional<
          Std::And( Capabilities::isAdaptiveDofMapper< Mapper >::v ... ),
          AdaptiveDofMapper< Traits< GridPart, Mapper ... > >,
          DofMapper< Traits< GridPart, Mapper ... > > >::type Type;
      };


    } // namespace __TupleMapper

#endif // #ifndef DOXYGEN


    // TupleMapper
    // -----------

    /** \brief mapper allocating one DoF per subentity of a given codimension
     *
     *  \tparam  GridPart  grid part type
     *  \tparam  ...Mapper  Parameter Pack of Mappers
     *
     *  \note This mapper is adaptve (cf. AdaptiveDofMapper) if and only if the
     *        grid part's index set is adaptive, i.e. if
     *        Capabilities::isAdaptiveIndexSet< GridPart::IndexSetType >::v is \b true
     */

    template< class GridPart, class ... Mapper >
    class TupleMapper
      : public __TupleMapper::template Implementation< GridPart, Mapper ... >::Type
    {
      typedef typename __TupleMapper::template Implementation< GridPart, Mapper ... >::Type BaseType;

    public:
      TupleMapper ( GridPart &gridPart, Mapper & ... mapper ) : BaseType( gridPart, mapper ... ) {}
      TupleMapper ( GridPart &gridPart, Mapper && ... mapper ) : BaseType( gridPart, std::move( mapper ) ... ) {}
    };

    // Capabilities
    // ------------

    namespace Capabilities
    {
      template< class GridPart, class ... Mapper >
      struct isAdaptiveDofMapper< TupleMapper< GridPart, Mapper ... > >
      {
        static const bool v = Std::And( isAdaptiveDofMapper< Mapper >::v ... );
      };

      template< class GridPart, class ... Mapper >
      struct isConsecutiveIndexSet< __TupleMapper::AdaptiveDofMapper< __TupleMapper::Traits< GridPart, Mapper ... > > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  }   // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH
