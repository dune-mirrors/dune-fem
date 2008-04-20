#ifndef DUNE_FEM_ATTACHEDFUNCTION_FUNCTION_HH
#define DUNE_FEM_ATTACHEDFUNCTION_FUNCTION_HH

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/function/common/dofblock.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>

#include "container.hh"

namespace Dune
{

  template< class DiscreteFunctionSpace >
  class AttachedDiscreteFunction;

  
  template< class DiscreteFunctionSpace >
  struct AttachedDiscreteFunctionTraits
  {
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType >
      DiscreteFunctionType;
    typedef AttachedDiscreteFunctionTraits< DiscreteFunctionSpaceType >
      Traits;

    typedef StandardLocalFunctionFactory< Traits > LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType >
      LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;
    
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType
      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

    typedef DofManager< typename GridPartType :: GridType > DofManagerType;

    typedef RangeFieldType DofType;
    typedef AttachedDiscreteFunctionContainer
      < DofType, DofManagerType, MapperType >
      ContainerType;

    typedef typename ContainerType :: SlotIteratorType DofIteratorType;
    typedef typename ContainerType :: ConstSlotIteratorType
      ConstDofIteratorType;

    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

    typedef DofBlockProxy< DiscreteFunctionType, DofType, blockSize >
      DofBlockType;
    typedef DofBlockProxy
      < const DiscreteFunctionType, const DofType, blockSize >
      ConstDofBlockType;
    typedef Envelope< DofBlockType > DofBlockPtrType;
    typedef Envelope< ConstDofBlockType > ConstDofBlockPtrType;
  };



  template< class DiscreteFunctionSpace >
  class AttachedDiscreteFunction
  : public DiscreteFunctionDefault
      < AttachedDiscreteFunctionTraits< DiscreteFunctionSpace > >
  {
    typedef AttachedDiscreteFunction< DiscreteFunctionSpace > ThisType;
    typedef DiscreteFunctionDefault
      < AttachedDiscreteFunctionTraits< DiscreteFunctionSpace > >
      BaseType;

  public:
    typedef AttachedDiscreteFunctionTraits< DiscreteFunctionSpace > Traits;

    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename Traits :: DomainFieldType DomainFieldType;
    typedef typename Traits :: RangeFieldType RangeFieldType;

    typedef typename Traits :: DomainType DomainType;
    typedef typename Traits :: RangeType RangeType;

    typedef typename Traits :: LocalFunctionFactoryType
      LocalFunctionFactoryType;

    typedef typename Traits :: DofManagerType DofManagerType;
    typedef typename Traits :: ContainerType ContainerType;

    typedef typename Traits :: DofIteratorType DofIteratorType;
    typedef typename Traits :: ConstDofIteratorType ConstDofIteratorType;

    typedef typename Traits :: DofBlockType DofBlockType;
    typedef typename Traits :: ConstDofBlockType ConstDofBlockType;
    typedef typename Traits :: DofBlockPtrType DofBlockPtrType;
    typedef typename Traits :: ConstDofBlockPtrType ConstDofBlockPtrType;

  protected:
    const LocalFunctionFactoryType lfFactory_;

    ContainerType &container_;
    unsigned int slot_;

  public:
    inline AttachedDiscreteFunction ( const std :: string &name,
                                      const DiscreteFunctionSpaceType &dfSpace )
    : BaseType( name, dfSpace, lfFactory_ ),
      lfFactory_( *this ),
      container_( ContainerType :: attach
        ( dofManager( dfSpace ), dfSpace.mapper() ) ),
      slot_( container_.allocSlot() )
    {}

    inline AttachedDiscreteFunction ( const ThisType &other )
    : BaseType( other.name(), other.space(), lfFactory_ ),
      lfFactory_( *this ),
      container_( ContainerType :: attach
        ( dofManager( other.space() ), other.space().mapper() ) ),
      slot_( container_.allocSlot() )
    {
      assign( other );
    }

    inline ~AttachedDiscreteFunction ()
    {
      container_.freeSlot( slot_ );
      ContainerType :: detach( container_ );
    }

    inline ConstDofIteratorType dbegin () const
    {
      return container().begin( slot_ );
    }

    inline DofIteratorType dbegin ()
    {
      return container().begin( slot_ );
    }

    inline ConstDofIteratorType dend () const
    {
      return container().end( slot_ );
    }

    inline DofIteratorType dend ()
    {
      return container().end( slot_ );
    }

    inline ConstDofBlockPtrType block ( unsigned int index ) const
    {
      typename ConstDofBlockType :: KeyType key( this, index );
      return ConstDofBlockPtrType( key );
    }

    inline DofBlockPtrType block ( unsigned int index )
    {
      typename DofBlockType :: KeyType key( this, index );
      return DofBlockPtrType( key );
    }

    inline const RangeFieldType &dof ( unsigned int index ) const
    {
      return container().dof( slot_, index );
    }

    inline RangeFieldType &dof ( unsigned int index )
    {
      return container().dof( slot_, index );
    }

    inline void enableDofCompression ()
    {
      container().enableDofCompression();
    }

    inline int size () const
    {
      return container_.size();
    }

    inline void swap ( ThisType &other )
    {
      const unsigned int myslot = slot_;
      slot_ = other.slot_;
      other.slot_ = myslot;
    }

  protected:
    inline const ContainerType &container () const
    {
      return container_;
    }

    inline ContainerType &container ()
    {
      return container_;
    }

  protected:
    inline static DofManagerType &
    dofManager ( const DiscreteFunctionSpaceType &dfSpace )
    {
      return
        DofManagerFactory< DofManagerType > :: getDofManager( dfSpace.grid() );
    }
  };

}

#endif
