#ifndef DUNE_FEM_NONBLOCKMAPPER_HH
#define DUNE_FEM_NONBLOCKMAPPER_HH

#include <dune/fem/space/common/dofmapper.hh>

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

    typedef typename BlockMapper::Traits::EntityType EntityType;

    typedef NonBlockDofMapIterator< typename BlockMapper::Traits::DofMapIteratorType, blockSize >
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

  public:
    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::DofMapIteratorType DofMapIteratorType;

    explicit NonBlockMapper ( BlockMapper &blockMapper )
    : blockMapper_( blockMapper )
    {}

    int size () const
    {
      return blockSize * blockMapper_.size();
    }

    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( blockMapper_.begin() );
    }

    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( blockMapper_.end() );
    }

    int mapToGlobal ( const EntityType &entity, const int localDof ) const
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

    int numDofs ( const EntityType &entity ) const
    {
      return blockSize * blockMapper_.numDofs( entity );
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

    int oldIndex ( const int hole, const int block )
    {
      const int i = hole % blockSize;
      const int blockHole = hole / blockSize;
      return blockMapper_.oldIndex( blockHole, block ) * blockSize + i;
    }

    int newIndex ( const int hole, const int block )
    {
      const int i = hole % blockSize;
      const int blockHole = hole / blockSize;
      return blockMapper_.newIndex( blockHole, block ) * blockSize + i;
    }

    bool consecutive () const
    {
      return blockMapper_.consecutive();
    }

    void update ( const bool oversize )
    {
      return blockMapper_.update( oversize );
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

  private:
    BlockMapper &blockMapper_;
  };

}

#endif // #ifndef DUNE_FEM_NONBLOCKMAPPER_HH
