#ifndef DUNE_FEM_NONBLOCKMAPPER_HH
#define DUNE_FEM_NONBLOCKMAPPER_HH

#include <vector>

#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class BlockMapper, int blockSize >
    class NonBlockMapper;


    namespace __NonBlockMapper
    {

      // Traits
      // ------

      template< class BlockMapper, int bS >
      struct Traits
      {
        typedef NonBlockMapper< BlockMapper, bS > DofMapperType;

        typedef BlockMapper BlockMapperType;
        typedef typename BlockMapper::ElementType ElementType;
        typedef typename BlockMapper::SizeType SizeType;
        typedef typename BlockMapper::GlobalKeyType GlobalKeyType;

        static const int blockSize = bS;
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

        typedef typename BaseType::Traits Traits;

        typedef typename Traits::BlockMapperType BlockMapperType;
        static const int blockSize = Traits::blockSize;

      public:
        typedef typename Traits::ElementType ElementType;
        typedef typename Traits::SizeType SizeType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

      private:
        template< class Functor >
        struct BlockFunctor
        {
          explicit BlockFunctor ( Functor functor )
            : functor_( functor )
          {}

          template< class GlobalKey >
          void operator() ( int localBlock, const GlobalKey globalKey )
          {
            int localDof = blockSize*localBlock;
            SizeType globalDof = blockSize*globalKey;
            const int localEnd = localDof + blockSize;
            while( localDof != localEnd )
              functor_( localDof++, globalDof++ );
          }

        private:
          Functor functor_;
        };

      public:

        DofMapper ( BlockMapperType &blockMapper )
          : blockMapper_( blockMapper )
        {}

        SizeType size () const { return blockSize * blockMapper_.size(); }

        bool contains ( const int codim ) const { return blockMapper_.contains( codim ); }

        bool fixedDataSize ( int codim ) const { return blockMapper_.fixedDataSize( codim ); }

        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const
        {
          blockMapper_.mapEach( element, BlockFunctor< Functor >( f ) );
        }

        void map ( const ElementType &element, std::vector< std::size_t > &indices ) const
        DUNE_DEPRECATED_MSG( "This function is out dated and will be deleted in version 1.5" )
        {
          indices.resize( numDofs( element ) );
          mapEach( element, AssignFunctor< std::vector< std::size_t > >( indices ) );
        }

        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const
        {
          blockMapper_.mapEachEntityDof( entity, BlockFunctor< Functor >( f ) );
        }

        int maxNumDofs () const { return blockSize * blockMapper_.maxNumDofs(); }

        SizeType numDofs ( const ElementType &element ) const { return blockSize * blockMapper_.numDofs( element ); }

        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const { return blockSize * blockMapper_.numEntityDofs( entity ); }

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

      private:
        BlockMapperType &blockMapper_;
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

        typedef typename BaseType::Traits Traits;

        typedef typename Traits::BlockMapperType BlockMapperType;
        static const int blockSize = Traits::blockSize;

        using BaseType::blockMapper_;

      public:
        typedef typename Traits::ElementType ElementType;
        typedef typename Traits::SizeType SizeType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

        AdaptiveDofMapper ( BlockMapperType &blockMapper )
          : BaseType( blockMapper )
        {}

        bool consecutive () const { return blockMapper_.consecutive(); }

        SizeType numBlocks () const { return blockMapper_.numBlocks(); }

        SizeType numberOfHoles ( const int block ) const { return blockSize * blockMapper_.numberOfHoles( block ); }

        GlobalKeyType oldIndex ( const int hole, const int block ) const
        {
          const int i = hole % blockSize;
          const int blockHole = hole / blockSize;
          return blockMapper_.oldIndex( blockHole, block ) * blockSize + i;
        }

        GlobalKeyType newIndex ( const int hole, const int block ) const
        {
          const int i = hole % blockSize;
          const int blockHole = hole / blockSize;
          return blockMapper_.newIndex( blockHole, block ) * blockSize + i;
        }

        SizeType oldOffSet ( const int block ) const { return blockMapper_.oldOffSet( block ) * blockSize; }

        SizeType offSet ( const int block ) const { return blockMapper_.offSet( block ) * blockSize; }
      };


      // Implementation
      // --------------

      template< class BlockMapper, int blockSize, bool adaptive = Capabilities::isAdaptiveDofMapper< BlockMapper >::v >
      class Implementation
      {
        typedef __NonBlockMapper::Traits< BlockMapper, blockSize > Traits;

      public:
        typedef typename std::conditional< adaptive,
                                           AdaptiveDofMapper< Traits >,
                                           DofMapper< Traits > >::type Type;
      };

    } // namespace __NonBlockMapper


    // NonBlockMapper
    // --------------

    /**Flatten the index-space of a given BlockMapper. */
    template< class BlockMapper, int blockSize >
    class NonBlockMapper
      : public __NonBlockMapper::template Implementation< BlockMapper, blockSize >::Type
    {
      typedef typename __NonBlockMapper::template Implementation< BlockMapper, blockSize >::Type BaseType;

    public:

      NonBlockMapper ( BlockMapper &blockMapper )
        : BaseType( blockMapper )
      {}
    };


    // NonBlockMapper for NonBlockMapper
    // ---------------------------------

    template< class BlockMapper, int innerBlockSize, int outerBlockSize >
    class NonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize >, outerBlockSize >
      : public NonBlockMapper< BlockMapper, innerBlockSize *outerBlockSize >
    {
      typedef NonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize >, outerBlockSize > ThisType;
      typedef NonBlockMapper< BlockMapper, innerBlockSize *outerBlockSize > BaseType;

    public:
      explicit NonBlockMapper ( const NonBlockMapper< BlockMapper, innerBlockSize > &blockMapper )
        : BaseType( blockMapper.blockMapper_ )
      {}
    };


    // Capabilities
    // ------------

    namespace Capabilities
    {
      template< class BlockMapper, int blockSize >
      struct isAdaptiveDofMapper< NonBlockMapper< BlockMapper, blockSize > >
      {
        static const bool v = isAdaptiveDofMapper< BlockMapper >::v;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NONBLOCKMAPPER_HH
