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

    // from adaptationmanager.hh
    template <class GridType>
    class AdaptationMethod ;


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

      protected:
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
          computeOffset();
        }

        DofMapper ( GridPartType &gridPart, Mapper && ... mapper )
          : gridPart_( gridPart ),
            mapperTuple_( std::move( mapper ) ... )
        {
          computeOffset();
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

          computeOffset();
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

        void computeOffset ()
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
        typedef typename GridPart :: GridType GridType;
        typedef AdaptationMethod< GridType > AdaptationMethodType;

        typedef typename BaseType::OffsetType OffsetType;

        using BaseType::mapperTupleSize;
        using BaseType::mapperTuple_;
        using BaseType::gridPart_;
        using BaseType::globalOffset_;

      public:
        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;
        typedef typename BaseType::GlobalKeyType GlobalKeyType;
        typedef GridPart GridPartType;

        AdaptiveDofMapper ( GridPartType &gridPart, Mapper & ... mapper )
          : BaseType( gridPart, mapper ... ),
            numBlocks_( computeNumBlocks() ),
            isCallBackAdapt_( AdaptationMethodType( gridPart.grid() ).isCallBackAdaptation() ),
            needsFullUpdate_(true)
        {
          oldGlobalOffset_ = globalOffset_;
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        AdaptiveDofMapper ( GridPartType &gridPart, Mapper && ... mapper )
          : BaseType( gridPart, std::move( mapper ) ... ),
            numBlocks_( computeNumBlocks() ),
            isCallBackAdapt_( AdaptationMethodType( gridPart.grid() ).isCallBackAdaptation() ),
            needsFullUpdate_(true)
        {
          oldGlobalOffset_ = globalOffset_;
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        ~AdaptiveDofMapper () { DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).removeIndexSet( *this ); }

        AdaptiveDofMapper ( const AdaptiveDofMapper & ) = delete;
        AdaptiveDofMapper ( AdaptiveDofMapper && ) = delete;

        static constexpr bool consecutive () noexcept { return true; }

        SizeType numBlocks () const { return numBlocks( std::index_sequence_for< Mapper ... >() ); }

        SizeType numberOfHoles ( const int block ) const
        {
          SizeType nHoles = 0;
            Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              const int localBlock = computeBlock( i, block );
              if( localBlock >= 0 )
              {
                nHoles = std::get< i >( this->mapperTuple_ ).numberOfHoles( localBlock );
                return;
              }
            } );
          return nHoles;
        }

        GlobalKeyType oldIndex ( const int hole, const int block ) const
        {
          SizeType oIndex = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              const int localBlock = computeBlock( i, block );
              if( localBlock >= 0 )
              {
                oIndex = std::get< i >( this->mapperTuple_ ).oldIndex( hole, localBlock ) + globalOffset_[ i ];
                return;
              }
            } );
          return oIndex;
        }

        GlobalKeyType newIndex ( const int hole, const int block ) const
        {
          SizeType nIndex = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              const int localBlock = computeBlock( i, block );
              if( localBlock >= 0 )
              {
                nIndex = std::get< i >( this->mapperTuple_ ).newIndex( hole, localBlock ) + globalOffset_[ i ];
                return ;
              }
            } );
          return nIndex;
        }

        SizeType oldOffSet ( int block ) const
        {
          SizeType oOffset = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              const int localBlock = computeBlock( i, block );
              if( localBlock >= 0 )
              {
                oOffset = std::get< i >( this->mapperTuple_ ).oldOffSet( localBlock ) + oldGlobalOffset_[ i ];
                return ;
              }
            } );
          return oOffset;
        }

        SizeType offSet ( int block ) const
        {
          SizeType offset = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i )
            {
              const int localBlock = computeBlock( i, block );
              if( localBlock >= 0 )
              {
                offset = std::get< i >( this->mapperTuple_ ).offSet( localBlock ) + globalOffset_[ i ];
                return ;
              }
            } );
          return offset;
        }

        void resize () { update(); }

        bool compress ()
        {
          update();
          // after this step for both methods we need the full update
          needsFullUpdate_ = true;
          return true;
        }

        void backup () const {}

        void restore () { compress(); }

        template< class IOStream >
        void read ( IOStream &in ) { compress(); }

        template< class IOStream >
        void write ( IOStream &out ) const {}

        template< class Entity >
        void insertEntity ( const Entity & )
        {
          if( needsFullUpdate_ )
          {
            // update also offsets in this case
            update();
          }
          else
          {
            // call BaseType::update to avoid changing oldGlobalOffSet_
            // do not use this->update here!
            BaseType::update();
          }
        }

        template< class Entity >
        void removeEntity ( const Entity & ) {}

        void update ()
        {
          // store previous offset
          oldGlobalOffset_ = globalOffset_;

          // in callback we always need the full update, otherwise set to false here
          needsFullUpdate_ = isCallBackAdapt_ ;

          // update component mappers and compute offsets
          BaseType::update();
        }

      protected:
        int computeBlock( const int i, const unsigned int block ) const
        {
          if( block >= blocks_[ i ] && block < blocks_[ i + 1 ] )
            return block - blocks_[ i ];
          else
            return -1;
        }

        void printOffSet( const OffsetType& offset, const std::string& name ) const
        {
          for( const auto& off : offset )
          {
            std::cout << name << " = " << off << std::endl;
          }
        }

        template< std::size_t ... i >
        SizeType numBlocks ( std::index_sequence< i ... > ) const
        {
          return Std::sum( std::get< i >( mapperTuple_ ).numBlocks() ... );
        }

        SizeType computeNumBlocks ()
        {
          // compute blocks (only needs to be done once)
          blocks_[ 0 ] = 0;
          Hybrid::forEach( std::make_index_sequence< mapperTupleSize >{},
            [ & ]( auto i ){ blocks_[ i + 1 ] = blocks_[ i ] + std::get< i >( mapperTuple_ ).numBlocks(); } );

          return blocks_[ mapperTupleSize + 1 ];
        }


        OffsetType oldGlobalOffset_;
        OffsetType blocks_;

        const SizeType numBlocks_;
        const bool isCallBackAdapt_;
        bool needsFullUpdate_;
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
