#ifndef DUNE_FEM_NONBLOCKMAPPER_HH
#define DUNE_FEM_NONBLOCKMAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class BlockMapper, int blockSize >
  class NonBlockMapper;



  // NonBlockDofMapIterator
  // ----------------------

  template< class BlockDofIterator, int blockSize >
  class NonBlockDofMapIterator
  {
    typedef NonBlockDofMapIterator< BlockDofIterator, blockSize > ThisType;

  public:
    NonBlockDofMapIterator ( const BlockDofIterator &blockIterator )
    : blockIterator_( blockIterator ), i_( 0 )
    {}

    ThisType &operator++ ()
    {
      if( ++i_ == blockSize )
      {
        i_ = 0;
        ++blockIterator_;
      }
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return (blockIterator_ == other.blockIterator_) && (i_ == other.i_);
    }

    bool operator!= ( const ThisType &other ) const
    {
      return (blockIterator_ != other.blockIterator_) || (i_ != other.i_);
    }

    int local () const
    {
      return blockSize * blockIterator_.local() + i_;
    }

    int global () const
    {
      return blockSize * blockIterator_.global() + i_;
    }

  private:
    BlockDofIterator blockIterator_;
    int i_;
  };



  // NonBlockMapperTraits
  // --------------------

  template< class BlockMapper, int blockSize >
  struct NonBlockMapperTraits
  {
    typedef NonBlockMapper< BlockMapper, blockSize > DofMapperType;

    typedef typename BlockMapper::ElementType ElementType;

    typedef NonBlockDofMapIterator< typename BlockMapper::DofMapIteratorType, blockSize >
      DofMapIteratorType;
  };



  // NonBlockMapper
  // --------------

  template< class BlockMapper, int blockSize >
  class NonBlockMapper
  : public DofMapper< NonBlockMapperTraits< BlockMapper, blockSize > >
  {
    typedef NonBlockMapper< BlockMapper, blockSize > ThisType;
    typedef DofMapper< NonBlockMapperTraits< BlockMapper, blockSize > > BaseType;

    template< class Functor >
    struct BlockFunctor
    {
      explicit BlockFunctor ( Functor functor )
      : functor_( functor )
      {}

      void operator() ( int localBlock, int globalBlock )
      {
        int localDof = blockSize*localBlock;
        int globalDof = blockSize*globalBlock;
        const int localEnd = localDof + blockSize;
        while( localDof != localEnd )
          functor_( localDof++, globalDof++ );
      }

    private:
      Functor functor_;
    };

  public:
    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::DofMapIteratorType DofMapIteratorType;

    explicit NonBlockMapper ( BlockMapper &blockMapper )
    : blockMapper_( blockMapper )
    {}

    int size () const
    {
      return blockSize * blockMapper_.size();
    }

    DofMapIteratorType begin ( const ElementType &entity ) const
    {
      return DofMapIteratorType( blockMapper_.begin(entity) );
    }

    DofMapIteratorType end ( const ElementType &entity ) const
    {
      return DofMapIteratorType( blockMapper_.end(entity) );
    }

    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const
    {
      blockMapper_.mapEach( element, BlockFunctor< Functor >( f ) );
    }

    int mapToGlobal ( const ElementType &entity, const int localDof ) const
    {
      const int i = localDof % blockSize;
      const int blockDof = localDof / blockSize;
      return blockMapper_.mapToGlobal( entity, blockDof ) * blockSize + i;
    }

    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      const int i = localDof % blockSize;
      const int blockDof = localDof / blockSize;
      return blockMapper_.mapEntityDofToGlobal( entity, blockDof ) * blockSize + i;
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

} // namespace Dune

#endif // #ifndef DUNE_FEM_NONBLOCKMAPPER_HH
