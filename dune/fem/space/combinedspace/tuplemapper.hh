#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH

#include <array>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>

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
          void operator() ( int localDof, const GlobalKey &globalKey ) const
          {
            functor_( localDof + localOffset_, CombinedIndex< GlobalKey, int, i >( globalKey, globalOffset_ ) );
          }

          template< class GlobalKey >
          void operator() ( const GlobalKey &globalKey ) const
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

        SizeType size () const { return size( std::index_sequence_for< Mapper ... >() ); }

        bool contains ( const int codim ) const { return contains( codim, std::index_sequence_for< Mapper ... >() ); }

        bool fixedDataSize ( int codim ) const { return fixedDataSize( codim, std::index_sequence_for< Mapper ... >() ); }

        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const
        {
          OffsetType localOffset;
          localOffset[ 0 ] = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              FunctorWrapper< Functor, i > wrappedFunctor( f, localOffset[ i ], globalOffset_[ i ] );
              std::get< i >( mapperTuple_ ).mapEach( element, wrappedFunctor );
              localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( mapperTuple_ ).numDofs( element );
            } );
        }

        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const
        {
          OffsetType localOffset;
          localOffset[ 0 ] = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              FunctorWrapper< Functor, i > wrappedFunctor( f, localOffset[ i ], globalOffset_[ i ] );
              std::get< i >( mapperTuple_ ).mapEachEntityDof( entity, wrappedFunctor );
              localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( mapperTuple_ ).numEntityDofs( entity );
            } );
        }

        void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
        {
          DUNE_THROW( NotImplemented, "Method onSubEntity(...) not yet implemented for TupleMapper" );
        }

        int maxNumDofs () const { return maxNumDofs( std::index_sequence_for< Mapper ... >() ); }

        SizeType numDofs ( const ElementType &element ) const { return numDofs( element, std::index_sequence_for< Mapper ... >() ); }

        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const { return numEntityDofs( entity, std::index_sequence_for< Mapper ... >() ); }


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

        void update()
        {
          // compute update for each mapper (if any)
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ](auto i){ std::get< i >( mapperTuple_ ).update(); } );
        }

        /*** NonInterface Methods ***/

        SizeType offset ( int i ) const { return globalOffset_[ i ]; }

        template< int i >
        SizeType subSize () const { return std::get< i >( mapperTuple_ ).size(); }

      protected:
        template< std::size_t ... i >
        SizeType size ( std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).size() ... );
        }

        template< std::size_t ... i >
        bool fixedDataSize ( const int codim, std::index_sequence< i ... > ) const
        {
          return Std::And( std::get< i >( mapperTuple_ ).fixedDataSize( codim ) ... );
        }

        template< std::size_t ... i >
        bool contains ( const int codim, std::index_sequence< i ... > ) const
        {
          return Std::Or( std::get< i >( mapperTuple_ ).contains( codim ) ... );
        }

        template< std::size_t ... i >
        int maxNumDofs ( std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).maxNumDofs() ... );
        }

        template< std::size_t ... i >
        SizeType numDofs ( const ElementType &element, std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numDofs( element ) ... );
        }

        template< class Entity, std::size_t ... i >
        SizeType numEntityDofs ( const Entity &entity, std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numEntityDofs( entity ) ... );
        }

        void init ()
        {
          globalOffset_[ 0 ] = 0;
          // compute new offsets
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i ){ globalOffset_[ i + 1 ] = globalOffset_[ i ] + std::get< i >( mapperTuple_ ).size(); } );
        }

        GridPartType &gridPart_;
        std::tuple< Mapper ... > mapperTuple_;
        OffsetType globalOffset_;
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

        SizeType numBlocks () const { return numBlocks( std::index_sequence_for< Mapper ... >() ); }

        SizeType numberOfHoles ( int block ) const
        {
          SizeType nHoles = 0;
          int comp = -1;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              if( comp >= 0 )
                return;
              const int localBlock = block - std::get< i >( this->mapperTuple_ ).numBlocks();
              if( localBlock >= 0 )
              {
                comp = i;
                nHoles = std::get< i >( this->mapperTuple_ ).numberOfHoles( localBlock );
              }
            } );
          return nHoles;
        }

        GlobalKeyType oldIndex ( int hole, int block ) const
        {
          int comp = -1;
          SizeType oIndex = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              if( comp >= 0 )
                return;
              const int localBlock = block - std::get< i >( this->mapperTuple_ ).numBlocks();
              if( localBlock >= 0 )
              {
                comp = i;
                oIndex = std::get< i >( this->mapperTuple_ ).oldIndex( hole, localBlock );
              }
            } );
          assert( comp >= 0 );
          return oIndex + globalOffset_[ comp ];
        }

        GlobalKeyType newIndex ( int hole, int block ) const
        {
          int comp = -1;
          SizeType nIndex = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              if( comp >= 0 )
                return;
              const int localBlock = block - std::get< i >( this->mapperTuple_ ).numBlocks();
              if( localBlock >= 0 )
              {
                comp = i;
                nIndex = std::get< i >( this->mapperTuple_ ).newIndex( hole, localBlock );
              }
            } );
          assert( comp > 0 );
          return nIndex + globalOffset_[ comp ];
        }

        SizeType oldOffSet ( int block ) const
        {
          int comp = -1;
          SizeType oOffset = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              if( comp >= 0 )
                return;
              const int localBlock = block - std::get< i >( this->mapperTuple_ ).numBlocks();
              if( localBlock >= 0 )
              {
                comp = i;
                oOffset = std::get< i >( this->mapperTuple_ ).oldOffSet( localBlock );
              }
            } );
          assert( comp >= 0 );
          return oOffset + oldGlobalOffset_[ comp ];
        }

        SizeType offSet ( int block ) const
        {
          int comp = -1;
          SizeType offset = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              if( comp >= 0 )
                return;
              const int localBlock = block - std::get< i >( this->mapperTuple_ ).numBlocks();
              if( localBlock >= 0 )
              {
                comp = i;
                offset = std::get< i >( this->mapperTuple_ ).offSet( localBlock );
              }
            } );
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

        void update ()
        {
          oldGlobalOffset_ = globalOffset_;
          BaseType::init();
        }

      protected:
        template< std::size_t ... i >
        SizeType numBlocks ( std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numBlocks() ... );
        }

      private:
        OffsetType oldGlobalOffset_;
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
