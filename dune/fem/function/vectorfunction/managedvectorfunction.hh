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
              class DofVector >
    class ManagedDiscreteFunction
      < VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
    : public VectorDiscreteFunction< DiscreteFunctionSpace, DofVector >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > BaseType;
      typedef ManagedDiscreteFunction< BaseType > ThisType;

    public:
      typedef ThisType DiscreteFunctionType;

      typedef typename BaseType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType DofVectorType;

    protected:
      typedef typename DiscreteFunctionSpaceType :: GridPartType :: GridType
        GridType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, 
                              DiscreteFunctionSpaceType::localBlockSize > NonBlockingMapperType ;

    protected:
      NonBlockingMapperType* mapper_ ;
      DofStorageInterface *memObject_;

    public:
      inline ManagedDiscreteFunction ( const std :: string &name,
                                       const DiscreteFunctionSpaceType &dfSpace )
      : BaseType( name, dfSpace, allocDofVector( name, dfSpace ) )
      {}

      inline explicit ManagedDiscreteFunction ( const BaseType &other )
      : BaseType( other.name(), other.space(),
                  allocDofVector( other.name(), other.space() ) )
      {
        BaseType :: assign ( other );
      }

      inline ManagedDiscreteFunction ( const ThisType &other )
      : BaseType( other.name(), other.space(),
                  allocDofVector( other.name(), other.space() ) )
      {
        BaseType :: assign ( other );
      }

      inline ~ManagedDiscreteFunction ()
      {
        if( memObject_ ) 
          delete memObject_ ;
        memObject_ = 0;

        delete mapper_ ; mapper_ = 0;
      }

      inline void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    private:
      inline DofVectorType &
      allocDofVector ( const std :: string &name,
                       const DiscreteFunctionSpaceType &dfSpace )
      {
        mapper_ = new NonBlockingMapperType( dfSpace.blockMapper() );
        // allocate managed dof storage 
        std::pair< DofStorageInterface *, DofVectorType * > memPair
          = allocateManagedDofStorage( dfSpace.gridPart().grid(), *mapper_, name, (DofVectorType *)0 );
        memObject_ = memPair.first;
        return *(memPair.second);
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
