#ifndef DUNE_FEM_NONBLOCKMAPPER_HH
#define DUNE_FEM_NONBLOCKMAPPER_HH

#include <vector>

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/misc/functor.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class BlockMapper, int blockSize >
    class NonBlockMapper;



    // NonBlockMapperTraits
    // --------------------

    template< class BlockMapper, int blockSize >
    struct NonBlockMapperTraits
    {
      typedef NonBlockMapper< BlockMapper, blockSize > DofMapperType;

      typedef typename BlockMapper::ElementType ElementType;

      typedef typename BlockMapper::SizeType  SizeType;
    };



    // NonBlockMapper
    // --------------

    /**Flatten the index-space of a given BlockMapper. */
    template< class BlockMapper, int blockSize >
    class NonBlockMapper
    : public DofMapper< NonBlockMapperTraits< BlockMapper, blockSize > >
    {
      typedef NonBlockMapper< BlockMapper, blockSize > ThisType;
      typedef DofMapper< NonBlockMapperTraits< BlockMapper, blockSize > > BaseType;

      template< class, int >
      friend class NonBlockMapper;

      typedef NonBlockMapperTraits< BlockMapper, blockSize > Traits ;
    public:
      typedef typename Traits :: SizeType SizeType ;

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
      typedef typename BaseType::ElementType ElementType;

      explicit NonBlockMapper ( BlockMapper &blockMapper )
      : blockMapper_( blockMapper )
      {}

      int size () const
      {
        return blockSize * blockMapper_.size();
      }

      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        blockMapper_.mapEach( element, BlockFunctor< Functor >( f ) );
      }

      void map ( const ElementType &element, std::vector< std::size_t > &indices ) const
      {
        indices.resize( numDofs( element ) );
        mapEach( element, AssignFunctor< std::vector< std::size_t > >( indices ) );
      }

      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        blockMapper_.mapEachEntityDof( entity, BlockFunctor< Functor >( f ) );
      }

      int maxNumDofs () const
      {
        return blockSize * blockMapper_.maxNumDofs();
      }

      int numDofs ( const ElementType &element ) const
      {
        return blockSize * blockMapper_.numDofs( element );
      }

      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        return blockSize * blockMapper_.numEntityDofs( entity );
      }

      int numberOfHoles ( const int block ) const
      {
        return blockSize * blockMapper_.numberOfHoles( block );
      }

      int oldIndex ( const int hole, const int block ) const
      {
        const int i = hole % blockSize;
        const int blockHole = hole / blockSize;
        return blockMapper_.oldIndex( blockHole, block ) * blockSize + i;
      }

      int newIndex ( const int hole, const int block ) const
      {
        const int i = hole % blockSize;
        const int blockHole = hole / blockSize;
        return blockMapper_.newIndex( blockHole, block ) * blockSize + i;
      }

      bool consecutive () const
      {
        return blockMapper_.consecutive();
      }

      int oldOffSet ( const int block ) const
      {
        return blockMapper_.oldOffSet( block ) * blockSize;
      }

      int offSet ( const int block ) const
      {
        return blockMapper_.offSet( block ) * blockSize;
      }

      int numBlocks () const
      {
        return blockMapper_.numBlocks();
      }

      bool contains( const int codim ) const
      {
        return blockMapper_.contains( codim );
      }

    private:
      BlockMapper &blockMapper_;
    };



    // NonBlockMapper for NonBlockMapper
    // ---------------------------------

    template< class BlockMapper, int innerBlockSize, int outerBlockSize >
    class NonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize >, outerBlockSize >
      : public NonBlockMapper< BlockMapper, innerBlockSize * outerBlockSize >
    {
      typedef NonBlockMapper< NonBlockMapper< BlockMapper, innerBlockSize >, outerBlockSize > ThisType;
      typedef NonBlockMapper< BlockMapper, innerBlockSize * outerBlockSize > BaseType;

    public:
      explicit NonBlockMapper ( const NonBlockMapper< BlockMapper, innerBlockSize > &blockMapper )
        : BaseType( blockMapper.blockMapper_ )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NONBLOCKMAPPER_HH
