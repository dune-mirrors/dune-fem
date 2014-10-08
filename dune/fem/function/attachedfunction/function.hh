#ifndef DUNE_FEM_ATTACHEDFUNCTION_FUNCTION_HH
#define DUNE_FEM_ATTACHEDFUNCTION_FUNCTION_HH

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/dofblock.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include "container.hh"

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunctionSpace >
    class AttachedDiscreteFunction;


    template< class DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< AttachedDiscreteFunction< DiscreteFunctionSpace > > 
    {
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType >
        DiscreteFunctionType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType
        DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType
        JacobianRangeType;

      enum { blockSize = DiscreteFunctionSpaceType::localBlockSize };

      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, blockSize > MapperType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef typename GridPartType::GridType GridType;

      typedef RangeFieldType DofType;
      typedef AttachedDiscreteFunctionContainer< DofType, GridType, MapperType >
        ContainerType;

      typedef typename ContainerType::SlotIteratorType DofIteratorType;
      typedef typename ContainerType::ConstSlotIteratorType
        ConstDofIteratorType;


      typedef DofBlockProxy< DiscreteFunctionType, DofType, blockSize >
        DofBlockType;
      typedef DofBlockProxy
        < const DiscreteFunctionType, const DofType, blockSize >
        ConstDofBlockType;
      typedef Fem::Envelope< DofBlockType > DofBlockPtrType;
      typedef Fem::Envelope< ConstDofBlockType > ConstDofBlockPtrType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    template< class DiscreteFunctionSpace >
    class AttachedDiscreteFunction
      : public DiscreteFunctionDefault< AttachedDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef AttachedDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< AttachedDiscreteFunction< DiscreteFunctionSpace > > BaseType;

    public:
      typedef DiscreteFunctionTraits< ThisType > Traits;

      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename BaseType::DomainFieldType DomainFieldType;
      typedef typename BaseType::RangeFieldType RangeFieldType;

      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;

      typedef typename Traits::ContainerType ContainerType;

      typedef typename BaseType::DofIteratorType DofIteratorType;
      typedef typename BaseType::ConstDofIteratorType ConstDofIteratorType;

      typedef typename BaseType::DofType DofType;
      typedef typename BaseType::DofBlockType DofBlockType;
      typedef typename BaseType::ConstDofBlockType ConstDofBlockType;
      typedef typename BaseType::DofBlockPtrType DofBlockPtrType;
      typedef typename BaseType::ConstDofBlockPtrType ConstDofBlockPtrType;

      typedef typename BaseType::LocalDofVectorAllocatorType LocalDofVectorAllocatorType;

      using BaseType::space;

    protected:
      typedef typename Traits::MapperType MapperType;

      typename Traits :: LocalDofVectorStackType ldvStack_;
      MapperType mapper_;

      ContainerType &container_;
      unsigned int slot_;

    public:
      inline AttachedDiscreteFunction ( const std::string &name,
                                        const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( std::max( sizeof( DofType ), sizeof( DofType* ) ) * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
          mapper_( dfSpace.blockMapper() ),
          container_( ContainerType::attach( dfSpace.gridPart().grid(), mapper_ ) ),
          slot_( container_.allocSlot() )
      {}

      inline AttachedDiscreteFunction ( const ThisType &other )
        : BaseType( other.name(), other.space(), LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( other.ldvStack_ ),
          mapper_( space().blockMapper() ),
          container_( ContainerType::attach( other.space().gridPart().grid(), mapper_ ) ),
          slot_( container_.allocSlot() )
      {
        assign( other );
      }

      inline ~AttachedDiscreteFunction ()
      {
        container_.freeSlot( slot_ );
        ContainerType::detach( container_ );
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
        typename ConstDofBlockType::KeyType key( this, index );
        return ConstDofBlockPtrType( key );
      }

      inline DofBlockPtrType block ( unsigned int index )
      {
        typename DofBlockType::KeyType key( this, index );
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
    };

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ATTACHEDFUNCTION_FUNCTION_HH
