#ifndef DUNE_FEM_ATTACHEDFUNCTION_CONTAINER_HH
#define DUNE_FEM_ATTACHEDFUNCTION_CONTAINER_HH

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/function/attachedfunction/storage.hh>
#include <dune/fem/function/attachedfunction/slotiterator.hh>

namespace Dune
{

  template< class Dof, class DofManager, class Mapper >
  class AttachedDiscreteFunctionContainer
  {
    typedef AttachedDiscreteFunctionContainer< Dof, DofManager, Mapper >
      ThisType;

    class SingletonKey;

    friend class DefaultSingletonFactory< SingletonKey, ThisType >;
    //friend class SingletonList< SingletonKey, ThisType >;

    typedef SingletonList< SingletonKey, ThisType > SingletonProvider;

  public:
    typedef Dof DofType;
    typedef DofManager DofManagerType;
    typedef Mapper MapperType;

  private:
    typedef AttachedDiscreteFunctionStorage< DofType > StorageType;

  public:
    typedef SlotIterator< typename StorageType :: IteratorType >
      SlotIteratorType;
    typedef SlotIterator< typename StorageType :: ConstIteratorType >
      ConstSlotIteratorType;

  private:
    DofManagerType &dofManager_;

    MemObjectInterface *memObject_;
    StorageType *storage_;

  private:
    inline explicit
    AttachedDiscreteFunctionContainer ( const SingletonKey &key )
    : dofManager_( key.dofManager() )
    {
      std :: pair< MemObjectInterface *, StorageType * > memPair
        = dofManager_.addDofSet( (StorageType *)0, key.mapper(), name() );
      memObject_ = memPair.first;
      storage_ = memPair.second;
    }

    // prohibit copying
    AttachedDiscreteFunctionContainer ( const ThisType & );

    inline ~AttachedDiscreteFunctionContainer ()
    {
      assert( memObject_ );
      dofManager_.removeDofSet( *memObject_ );
    }

    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline unsigned int allocSlot ()
    {
      return storage().allocDof();
    }

    inline void freeSlot ( unsigned int slot )
    {
      return storage().freeDof( slot );
    }

    inline ConstSlotIteratorType begin ( unsigned int slot ) const
    {
      return ConstSlotIteratorType( storage().begin(), slot );
    }

    inline SlotIteratorType begin ( unsigned int slot )
    {
      return SlotIteratorType( storage().begin(), slot );
    }

    inline ConstSlotIteratorType end ( unsigned int slot ) const
    {
      return ConstSlotIteratorType( storage().end(), slot );
    }

    inline SlotIteratorType end ( unsigned int slot )
    {
      return SlotIteratorType( storage().end(), slot );
    }

    inline const DofType &dof ( unsigned int slot,
                                unsigned int index ) const
    {
      return storage()[ index ][ slot ];
    }

    inline DofType &dof ( unsigned int slot,
                          unsigned int index )
    {
      return storage()[ index ][ slot ];
    }

    inline void enableDofCompression ()
    {
      assert( memObject_ );
      memObject_->enableDofCompression();
    }

    inline const std :: string &name () const
    {
      static std :: string staticName( "AttachedDiscreteFunctionContainer" );
      return staticName;
    }

    inline unsigned int size () const
    {
      return storage().size();
    }

  protected:
    inline const StorageType &storage () const
    {
      assert( storage_ != 0 );
      return *storage_;
    }

    inline StorageType &storage ()
    {
      assert( storage_ != 0 );
      return *storage_;
    }

  public:
    inline static ThisType &attach ( DofManagerType &dofManager,
                                     MapperType &mapper )
    {
      SingletonKey key( dofManager, mapper );
      return SingletonProvider :: getObject( key );
    }

    inline static void detach ( const ThisType &instance )
    {
      SingletonProvider :: removeObject( instance );
    }
  };



  template< class Dof, class DofManager, class Mapper >
  class AttachedDiscreteFunctionContainer< Dof, DofManager, Mapper >
    :: SingletonKey
  {
    typedef SingletonKey ThisType;

  public:
    typedef DofManager DofManagerType;
    typedef Mapper MapperType;

  protected:
    DofManager *const dofManager_;
    MapperType *const mapper_;

  public:
    inline SingletonKey ( DofManagerType &dofManager,
                          MapperType &mapper )
    : dofManager_( &dofManager ),
      mapper_( &mapper )
    {}

    inline bool operator== ( const ThisType &other ) const
    {
      return ((dofManager_ == other.dofManager_) && (mapper_ == other.mapper_));
    }

    inline DofManagerType &dofManager () const
    {
      return *dofManager_;
    }

    inline MapperType &mapper () const
    {
      return *mapper_;
    }
  };

}

#endif
