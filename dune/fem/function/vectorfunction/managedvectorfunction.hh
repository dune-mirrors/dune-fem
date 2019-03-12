#ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
#define DUNE_FEM_MANAGEDVECTORFUNCTION_HH

#include <memory>
#include <string>

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // MemObject
    // ---------

    struct MemObject
    {
      template< class DofContainer, class DiscreteFunctionSpace >
      DofContainer &allocate ( const DiscreteFunctionSpace &space )
      {
        typedef MutableBlockVector< DofContainer, DiscreteFunctionSpace::localBlockSize > DofVector;
        auto result = allocateManagedDofStorage( space.gridPart().grid(), space.blockMapper(), static_cast< DofVector * >( nullptr ) );
        interface_.reset( result.first );
        return result.second->array();
      }

      void enableDofCompression ()
      {
        if( interface_ )
          interface_->enableDofCompression();
      }

    private:
      std::unique_ptr< DofStorageInterface > interface_;
    };




    // ManagedDiscreteFunction
    // -----------------------

    template< class DiscreteFunctionSpace,
              class Vector >
    class ManagedDiscreteFunction
      < VectorDiscreteFunction< DiscreteFunctionSpace, Vector > >
      : private MemObject,
        public VectorDiscreteFunction< DiscreteFunctionSpace, Vector >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector > BaseType;
      typedef ManagedDiscreteFunction< BaseType > ThisType;

    public:
      typedef ThisType DiscreteFunctionType;

      typedef typename BaseType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType    DofVectorType;
      typedef typename BaseType :: DofContainerType DofContainerType;

      typedef typename DiscreteFunctionSpaceType :: GridPartType :: GridType GridType;

      using BaseType :: name;

      ManagedDiscreteFunction ( const std::string &name, const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( name, dfSpace, MemObject::allocate< DofContainerType >( dfSpace ) )
      {}

      explicit ManagedDiscreteFunction ( const BaseType &other )
        : BaseType( other.name(), other.space(), MemObject::allocate< DofContainerType >( other.space() ) )
      {
        BaseType :: assign ( other );
      }

      ManagedDiscreteFunction ( const ThisType &other )
        : BaseType( other.name(), other.space(), MemObject::allocate< DofContainerType >( other.space() ) )
      {
        BaseType :: assign ( other );
      }

      ManagedDiscreteFunction ( ThisType && ) = default;

      void enableDofCompression () { MemObject::enableDofCompression(); }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
