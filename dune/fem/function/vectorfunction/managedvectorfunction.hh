#ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
#define DUNE_FEM_MANAGEDVECTORFUNCTION_HH

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunctionSpace,
              class Vector >
    class ManagedDiscreteFunction
      < VectorDiscreteFunction< DiscreteFunctionSpace, Vector > >
    : public VectorDiscreteFunction< DiscreteFunctionSpace, Vector >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector > BaseType;
      typedef ManagedDiscreteFunction< BaseType > ThisType;

    public:
      typedef ThisType DiscreteFunctionType;

      typedef typename BaseType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType    DofVectorType;
      typedef typename BaseType :: DofContainerType DofContainerType;

    protected:
      typedef typename DiscreteFunctionSpaceType :: GridPartType :: GridType
        GridType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef NonBlockMapper< BlockMapperType,
                              DiscreteFunctionSpaceType::localBlockSize > NonBlockingMapperType ;

    public:
      inline ManagedDiscreteFunction ( const std :: string &name,
                                       const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( name, dfSpace, allocDofContainer( name, dfSpace ) )
      {}

      inline explicit ManagedDiscreteFunction ( const BaseType &other )
      : BaseType( other.name(), other.space(),
                  allocDofContainer( other.name(), other.space() ) )
      {
        BaseType :: assign ( other );
      }

      inline ManagedDiscreteFunction ( const ThisType &other )
      : BaseType( other.name(), other.space(),
                  allocDofContainer( other.name(), other.space() ) )
      {
        BaseType :: assign ( other );
      }

      inline ~ManagedDiscreteFunction ()
      {
        if( memObject_ )
          delete memObject_ ;
        memObject_ = 0;
      }

      inline void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:
      inline DofContainerType &
      allocDofContainer ( const std :: string &name,
                          const DiscreteFunctionSpaceType &space )
      {
        typedef MutableBlockVector< DofContainerType,
                                    DiscreteFunctionSpaceType::localBlockSize > MutableDofVectorType;

        // allocate managed dof storage
        std::pair< DofStorageInterface *, MutableDofVectorType* > memPair
          = allocateManagedDofStorage( space.gridPart().grid(), space.blockMapper(), name, (MutableDofVectorType *)0 );
        memObject_ = memPair.first;
        return memPair.second->array();
      }

      // pointer to memory if allocated locally
      DofStorageInterface *memObject_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
