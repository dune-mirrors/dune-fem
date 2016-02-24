#ifndef DUNE_FEM_COMBINEDSPACE_POWERMAPPER_HH
#define DUNE_FEM_COMBINEDSPACE_POWERMAPPER_HH

#include <dune/common/math.hh>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {
    // Internal Forward Declration
    // ---------------------------

    template< class GridPart, class Mapper, int N >
    class PowerMapper;

#ifndef DOXYGEN

    namespace __PowerMapper
    {

      // Traits
      // ------

      template< class GridPart, class Mapper, int N >
      struct Traits
      {
        typedef Mapper MapperType;

        typedef GridPart GridPartType;

        typedef typename MapperType::ElementType ElementType;
        typedef typename MapperType::SizeType SizeType;
        typedef typename MapperType::GlobalKeyType GlobalKeyType;

        typedef PowerMapper< GridPartType, MapperType, N >  DofMapperType;

        static const int numComponents = N;
      };



      template< class T, template< class > class Base = Dune::Fem::DofMapper >
      class DofMapper
        : public Base< T >
      {
        typedef Base< T > BaseType;

      public:
        typedef typename BaseType::Traits Traits;
        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;

        typedef typename Traits::GridPartType GridPartType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

        typedef typename Traits::MapperType MapperType;

        static const int numComponents = Traits::numComponents;

      protected:
        // funtor wrapper
        template< class Functor >
        struct FunctorWrapper
        {
          FunctorWrapper ( SizeType offset, Functor functor )
            : offset_( offset ),
            functor_( functor )
          {}

          template< class GlobalKey >
          void operator() ( int localBlock, const GlobalKey &globalKey )
          {
            int localDof = localBlock*numComponents;
            for( int component = 0; component < numComponents; ++component, ++localDof )
              functor_( localDof, globalKey + component * offset_ );
          }

          template< class GlobalKey >
          void operator() ( const GlobalKey &globalKey )
          {
            for( int component = 0; component < numComponents; ++component )
              functor_( globalKey + component * offset_ );
          }

        private:
          SizeType offset_;
          Functor functor_;
        };

      public:

        DofMapper ( GridPartType &gridPart, MapperType &mapper )
          : gridPart_( gridPart ),
          mapper_( mapper ),
          offset_( mapper_.size() )
        {}

        SizeType size () const {  return mapper().size() * numComponents; }

        bool contains ( const int codim ) const { return mapper().contains( codim ); }

        bool fixedDataSize ( const int codim ) const { return mapper().fixedDataSize( codim ); }

        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const
        {
          FunctorWrapper< Functor > wrapper( offset_, f );
          mapper().mapEach( element, wrapper );
        }

        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const
        {
          FunctorWrapper< Functor > wrapper( offset_, f );
          mapper().mapEachEntityDof( entity, wrapper );
        }

        int maxNumDofs () const { return mapper().maxNumDofs() * numComponents; }

        SizeType numDofs ( const ElementType &element ) const { return mapper().numDofs( element ) * numComponents; }

        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const { return mapper().numEntityDofs( entity ) * numComponents; }


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

      protected:
        MapperType &mapper () { return mapper_; }
        const MapperType &mapper () const { return mapper_; }

        GridPartType &gridPart_;
        MapperType &mapper_;
        SizeType offset_;
      };



      // AdaptiveDofMapper
      // -----------------

      template< class T >
      class AdaptiveDofMapper
        : public DofMapper< T, Dune::Fem::AdaptiveDofMapper >
      {
        typedef DofMapper< T, Dune::Fem::AdaptiveDofMapper > BaseType;

      public:
        typedef typename BaseType::Traits Traits;
        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;

        typedef typename Traits::GridPartType GridPartType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

        typedef typename Traits::MapperType MapperType;

        static const int numComponents = Traits::numComponents;

      protected:
        using BaseType::mapper;
        using BaseType::gridPart_;
        using BaseType::offset_;

      public:

        AdaptiveDofMapper ( GridPartType &gridPart, MapperType &mapper )
          : BaseType( gridPart, mapper ),
          oldOffset_( -1 )
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        ~AdaptiveDofMapper ()
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).removeIndexSet( *this );
        }

        bool consecutive () const { return true; }

        SizeType numBlocks () const { return mapper().numBlocks() * numComponents; }

        SizeType numberOfHoles ( const int block ) const { return mapper().numberOfHoles( block % mapper().numBlocks() ); }

        GlobalKeyType oldIndex ( const int hole, const int block ) const
        {
          const int numContainedBlocks = mapper().numBlocks();
          const int containedBlock     = block % numContainedBlocks;
          const int component          = block / numContainedBlocks;

          const int containedOffset = mapper().oldIndex( hole, containedBlock );
          return containedOffset + component * oldOffset_;
        }

        GlobalKeyType newIndex ( const int hole, const int block ) const
        {
          const int numContainedBlocks = mapper().numBlocks();
          const int containedBlock     = block % numContainedBlocks;
          const int component          = block / numContainedBlocks;

          const int containedOffset = mapper().newIndex( hole, containedBlock );
          return containedOffset + component * offset_;
        }

        SizeType oldOffSet ( const int block ) const
        {
          const int numContainedBlocks = mapper().numBlocks();
          const int containedBlock     = block % numContainedBlocks;
          const int component          = block / numContainedBlocks;

          const int containedOffset = mapper().oldOffSet( containedBlock );
          return containedOffset + component * oldOffset_;
        }

        SizeType offSet ( const int block ) const
        {
          const int numContainedBlocks = mapper().numBlocks();
          const int containedBlock     = block % numContainedBlocks;
          const int component          = block / numContainedBlocks;

          const int containedOffset = mapper().offSet( containedBlock );
          return containedOffset + component * offset_;
        }

        template< class Entity >
        void insertEntity ( const Entity &entity ) { update(); }

        template< class Entity >
        void removeEntity ( const Entity &entity ) { }

        void resize () { update(); }

        bool compress () { update(); return true; }

        template< class StreamTraits >
        void write ( OutStreamInterface< StreamTraits > &out ) const {}

        template< class StreamTraits >
        void read ( InStreamInterface< StreamTraits > &in )
        {
          update();
        }

        void backup () const {}
        void restore () {}

      protected:
        void update ()
        {
          oldOffset_ = offset_;
          offset_ = mapper().size();
        }

      private:
        SizeType oldOffset_;
      };


      // Implementation
      // --------------

      template< class GridPart, class Mapper, int N, bool adapative = Capabilities::isAdaptiveDofMapper< Mapper >::v >
      class Implementation
      {
        typedef __PowerMapper::Traits< GridPart, Mapper, N > Traits;

      public:
        typedef typename std::conditional< adapative, AdaptiveDofMapper< Traits >, DofMapper< Traits > >::type Type;
      };



    } // namespace __PowerMapper

#endif // #ifndef DOXYGEN


    // PowerMapper
    // ----------

    /** \brief mapper allocating one DoF per subentity of a given codimension
     *
     *  \tparam  GridPart  grid part type
     *  \tparam  Mapper   contained mapper type
     *  \tparam  N  number of containd components
     *
     *  \note This mapper is adaptve (cf. AdaptiveDofMapper) if and only if the
     *        underlaying mapper is adaptiv, too.
     *        Capabilities::isAdaptiveDofMapper< Mapper >::v is \b true
     */

    template< class GridPart, class Mapper, int N >
    class PowerMapper
      : public __PowerMapper::template Implementation< GridPart, Mapper, N >::Type
    {
      typedef typename __PowerMapper::template Implementation< GridPart, Mapper, N >::Type BaseType;

    public:
      PowerMapper ( GridPart &gridPart, Mapper &mapper )
        : BaseType( gridPart, mapper )
      {}
    };


    // Capabilities
    // ------------

    namespace Capabilities
    {
      template< class GridPart, class Mapper, int N >
      struct isAdaptiveDofMapper< PowerMapper< GridPart, Mapper, N > >
      {
        static const bool v = isAdaptiveDofMapper< Mapper >::v;
      };

      template< class GridPart, class Mapper, int N >
      struct isConsecutiveIndexSet< __PowerMapper::AdaptiveDofMapper< __PowerMapper::Traits< GridPart, Mapper, N > > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities
      /** @} **/

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_COMBINEDSPACE_POWERMAPPER_HH
