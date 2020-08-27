#ifndef DUNE_FEM_FUNCTION_HIERARCHICAL_FUNCTION_HH
#define DUNE_FEM_FUNCTION_HIERARCHICAL_FUNCTION_HH

#include <dune/common/fvector.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/hierarchical/dofvector.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declaration
    // ----------------------------

    template< class DiscreteFunctionSpace >
    class HierarchicalDiscreteFunction;



    namespace Impl
    {

      template< class Dof, class BlockIndices >
      struct HierarchicalDofContainerChooser;

#if HAVE_DUNE_ISTL
      template< class Dof, int sz >
      struct HierarchicalDofContainerChooser< Dof, Hybrid::IndexRange< int, sz > >
      {
        typedef BlockVector< FieldVector< Dof, sz > > Type;
      };

      template< class Dof, class... SR >
      struct HierarchicalDofContainerChooser< Dof, Hybrid::CompositeIndexRange< SR... > >
      {
        typedef MultiTypeBlockVector< typename HierarchicalDofContainerChooser< Dof, SR >::Type... > Type;
      };
#else
      template< class Dof, int sz >
      struct HierarchicalDofContainerChooser< Dof, Hybrid::IndexRange< int, sz > >
      {
        typedef MutableBlockVector< DynamicArray< Dof >, sz > Type;
      };

      template< class Dof, class... SR >
      struct HierarchicalDofContainerChooser< Dof, Hybrid::CompositeIndexRange< SR... > >
      {
        typedef MutableBlockVector< DynamicArray< Dof >, Hybrid::CompositeIndexRange< SR... >::size() > Type;
      };
#endif // #if HAVE_DUNE_ISTL

    } // namespace Impl



    // DiscreteFunctionTraits for HierarchicalDiscreteFunction
    // -------------------------------------------------------

    template< class DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< HierarchicalDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::RangeFieldType DofType;

      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      typedef HierarchicalDofVector< typename Impl::HierarchicalDofContainerChooser< DofType, typename DiscreteFunctionSpaceType::LocalBlockIndices >::Type > DofVectorType;

      // fake DoF blocks, DoF block pointers and DoF iterators
      typedef DofType *DofIteratorType;
      typedef const DofType *ConstDofIteratorType;
      typedef DofType *DofBlockType;
      typedef const DofType *ConstDofBlockType;
      typedef DofType **DofBlockPtrType;
      typedef const DofType **ConstDofBlockPtrType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType * > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef HierarchicalDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    // HierarchicalDiscreteFunction
    // ----------------------------

    template< class DiscreteFunctionSpace >
    class HierarchicalDiscreteFunction
      : public DiscreteFunctionDefault< HierarchicalDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef HierarchicalDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< HierarchicalDiscreteFunction< DiscreteFunctionSpace > > BaseType;

    public:
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename BaseType::DofVectorType DofVectorType;

      using BaseType :: name;

      HierarchicalDiscreteFunction ( const std::string &name, const DiscreteFunctionSpaceType &space, DofVectorType &dofVector )
        : BaseType( name, space ), dofVector_( dofVector )
      {}

      HierarchicalDiscreteFunction ( const std::string &name, const DiscreteFunctionSpaceType &space )
        : BaseType( name, space ), dofVector_( allocateDofStorage( space ) )
      {}

      HierarchicalDiscreteFunction ( const ThisType &other )
        : BaseType( "copy of " + other.name(), other.space() ), dofVector_( allocateDofStorage( other.space() ) )
      {
        dofVector_ = other.dofVector_;
      }

      DofVectorType &dofVector () { return dofVector_; }
      const DofVectorType &dofVector () const { return dofVector_; }

      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:
      DofVectorType &allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        auto memPair = allocateManagedDofStorage< DofVectorType >( space.gridPart().grid(), space.blockMapper() );
        memObject_.reset( memPair.first );
        return *memPair.second;
      }

      std::unique_ptr< DofStorageInterface > memObject_;
      DofVectorType &dofVector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HIERARCHICAL_FUNCTION_HH
