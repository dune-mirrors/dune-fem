#ifndef DUNE_FEM_DYNAMICNONBLOCKMAPPER_HH
#define DUNE_FEM_DYNAMICNONBLOCKMAPPER_HH

#include <vector>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class BlockMapper >
    class DynamicNonBlockMapper;


    namespace __DynamicNonBlockMapper
    {

      // Traits
      // ------

      template< class BlockMapper >
      struct Traits
      {
        typedef DynamicNonBlockMapper< BlockMapper > DofMapperType;

        typedef BlockMapper BlockMapperType;
        typedef typename BlockMapper::ElementType ElementType;
        typedef typename BlockMapper::SizeType SizeType;
        typedef typename BlockMapper::GlobalKeyType GlobalKeyType;
      };


      // DofMapper
      // ---------

      template< class T, template< class > class Base = Dune::Fem::DofMapper >
      class DofMapper
        : public Base< T >
      {
        typedef Base< T > BaseType;

        template< class, template< class > class >
        friend class DofMapper;

      public:
        typedef typename BaseType::Traits Traits;

        typedef typename Traits::BlockMapperType BlockMapperType;

        typedef typename Traits::ElementType ElementType;
        typedef typename Traits::SizeType SizeType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

      private:
        template< class Functor >
        struct BlockFunctor
        {
          explicit BlockFunctor ( int blockSize, Functor functor )
            : blockSize_( blockSize ), functor_( functor )
          {}

          template< class GlobalKey >
          void operator() ( int localBlock, const GlobalKey globalKey )
          {
            int localDof = blockSize_*localBlock;
            SizeType globalDof = blockSize_*globalKey;
            const int localEnd = localDof + blockSize_;
            while( localDof != localEnd )
              functor_( localDof++, globalDof++ );
          }

        private:
          int blockSize_;
          Functor functor_;
        };

      public:
        DofMapper ( BlockMapperType &blockMapper, int blockSize )
          : blockMapper_( blockMapper ), blockSize_( blockSize )
        {}

        SizeType size () const { return blockSize() * blockMapper_.size(); }

        bool contains ( const int codim ) const { return blockMapper_.contains( codim ); }

        bool fixedDataSize ( int codim ) const { return blockMapper_.fixedDataSize( codim ); }

        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const
        {
          blockMapper_.mapEach( element, BlockFunctor< Functor >( blockSize(), f ) );
        }

        void map ( const ElementType &element, std::vector< GlobalKeyType > &indices ) const
        {
          indices.resize( numDofs( element ) );
          mapEach( element, [ &indices ] ( int local, GlobalKeyType global ) { indices[ local ] = global; } );
        }

        void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
        {
          const SizeType numDofs = blockMapper_.numDofs( element );
          blockMapper_.onSubEntity( element, i, c, indices );
          indices.resize( blockSize() * numDofs );
          for( SizeType i = numDofs; i > 0; )
          {
            for( int j = 0; j < blockSize(); ++j )
              indices[ i*blockSize() + j ] = indices[ i ];
          }
        }

        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const
        {
          blockMapper_.mapEachEntityDof( entity, BlockFunctor< Functor >( blockSize(), f ) );
        }

        template< class Entity >
        void mapEntityDofs ( const Entity &entity, std::vector< GlobalKeyType > &indices ) const
        {
          indices.resize( numEntityDofs( entity ) );
          mapEachEntityDof( entity, [ &indices ] ( int local, GlobalKeyType global ) { indices[ local ] = global; } );
        }

        int maxNumDofs () const { return blockSize() * blockMapper_.maxNumDofs(); }

        SizeType numDofs ( const ElementType &element ) const { return blockSize() * blockMapper_.numDofs( element ); }

        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const { return blockSize() * blockMapper_.numEntityDofs( entity ); }

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

        const BlockMapperType &blockMapper () const { return blockMapper_; }
        int blockSize () const { return blockSize_; }

      private:
        BlockMapperType &blockMapper_;
        int blockSize_;
      };


      // AdaptiveDofMapper
      // -----------------

      template< class T >
      class AdaptiveDofMapper
        : public DofMapper< T, Dune::Fem::AdaptiveDofMapper >
      {
        typedef DofMapper< T, Dune::Fem::AdaptiveDofMapper > BaseType;

        template< class >
        friend class AdaptiveDofMapper;

      public:
        typedef typename BaseType::Traits Traits;

        typedef typename Traits::BlockMapperType BlockMapperType;

        using BaseType::blockMapper;
        using BaseType::blockSize;

        typedef typename Traits::ElementType ElementType;
        typedef typename Traits::SizeType SizeType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

        AdaptiveDofMapper ( BlockMapperType &blockMapper, int blockSize )
          : BaseType( blockMapper, blockSize )
        {}

        bool consecutive () const { return blockMapper().consecutive(); }

        SizeType numBlocks () const { return blockMapper().numBlocks(); }

        SizeType numberOfHoles ( int block ) const { return blockSize() * blockMapper().numberOfHoles( block ); }

        GlobalKeyType oldIndex ( int hole, int block ) const
        {
          const int i = hole % blockSize();
          const int blockHole = hole / blockSize();
          return blockMapper().oldIndex( blockHole, block ) * blockSize() + i;
        }

        GlobalKeyType newIndex ( int hole, int block ) const
        {
          const int i = hole % blockSize;
          const int blockHole = hole / blockSize();
          return blockMapper().newIndex( blockHole, block ) * blockSize() + i;
        }

        SizeType oldOffSet ( const int block ) const { return blockMapper().oldOffSet( block ) * blockSize(); }

        SizeType offSet ( const int block ) const { return blockMapper().offSet( block ) * blockSize(); }
      };


      // Implementation
      // --------------

      template< class BlockMapper, bool adaptive = Capabilities::isAdaptiveDofMapper< BlockMapper >::v >
      class Implementation
      {
        typedef __DynamicNonBlockMapper::Traits< BlockMapper > Traits;

      public:
        typedef typename std::conditional< adaptive, AdaptiveDofMapper< Traits >, DofMapper< Traits > >::type Type;
      };

    } // namespace __DynamicNonBlockMapper



    // DynamicNonBlockMapper
    // ---------------------

    /** Flatten the index-space of a given BlockMapper **/
    template< class BlockMapper >
    class DynamicNonBlockMapper
      : public __DynamicNonBlockMapper::template Implementation< BlockMapper >::Type
    {
      typedef typename __DynamicNonBlockMapper::template Implementation< BlockMapper >::Type BaseType;

    public:
      DynamicNonBlockMapper ( BlockMapper &blockMapper, int blockSize )
        : BaseType( blockMapper, blockSize )
      {}
    };



    // DynamicNonBlockMapper for DynamicNonBlockMapper
    // -----------------------------------------------

    template< class BlockMapper >
    class DynamicNonBlockMapper< DynamicNonBlockMapper< BlockMapper > >
      : public DynamicNonBlockMapper< BlockMapper >
    {
      typedef DynamicNonBlockMapper< DynamicNonBlockMapper< BlockMapper > > ThisType;
      typedef DynamicNonBlockMapper< BlockMapper > BaseType;

    public:
      explicit DynamicNonBlockMapper ( const DynamicNonBlockMapper< BlockMapper > &blockMapper, int blockSize )
        : BaseType( blockMapper.blockMapper(), blockMapper.blockSize() * blockSize )
      {}
    };



    // DynamicNonBlockMapper for NonBlockMapper
    // ----------------------------------------

    template< class BlockMapper, int innerBlockSize >
    class DynamicNonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize > >
      : public DynamicNonBlockMapper< BlockMapper >
    {
      typedef DynamicNonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize > > ThisType;
      typedef DynamicNonBlockMapper< BlockMapper > BaseType;

    public:
      explicit DynamicNonBlockMapper ( const NonBlockMapper< BlockMapper, innerBlockSize > &blockMapper, int blockSize )
        : BaseType( blockMapper.blockMapper(), innerBlockSize * blockSize )
      {}
    };



    // NonBlockMapper for DynamicNonBlockMapper
    // ----------------------------------------

    template< class BlockMapper, int outerBlockSize >
    class NonBlockMapper< DynamicNonBlockMapper< BlockMapper >, outerBlockSize >
      : public DynamicNonBlockMapper< BlockMapper >
    {
      typedef NonBlockMapper< DynamicNonBlockMapper< BlockMapper >, outerBlockSize > ThisType;
      typedef DynamicNonBlockMapper< BlockMapper > BaseType;

    public:
      explicit NonBlockMapper ( const DynamicNonBlockMapper< BlockMapper > &blockMapper )
        : BaseType( blockMapper.blockMapper(), outerBlockSize * blockMapper.blockSize() )
      {}
    };



    // Capabilities
    // ------------

    namespace Capabilities
    {
      template< class BlockMapper >
      struct isAdaptiveDofMapper< DynamicNonBlockMapper< BlockMapper > >
      {
        static const bool v = isAdaptiveDofMapper< BlockMapper >::v;
      };

      template< class BlockMapper >
      struct isConsecutiveIndexSet< __DynamicNonBlockMapper::AdaptiveDofMapper< __DynamicNonBlockMapper::Traits< BlockMapper > > >
      {
        static const bool v = isConsecutiveIndexSet< BlockMapper >::v;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DYNAMICNONBLOCKMAPPER_HH
