#ifndef DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH
#define DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH

#include <dune/fem/function/common/scalarproducts.hh>


#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>


namespace Dune
{

  namespace Fem
  {

    /* =================================================
     * class PetscSlaveDofProvider
     * =================================================
     */
    template< typename DFSpace >
    class PetscSlaveDofProvider : public SlaveDofsProvider< DFSpace >
    {
    public:
      typedef DFSpace DiscreteFunctionSpaceType;

      typedef PetscSlaveDofProvider< DiscreteFunctionSpaceType > ThisType;
      typedef SlaveDofsProvider< DiscreteFunctionSpaceType >     BaseType;
    protected:
      using BaseType :: space_;
      using BaseType :: slaveDofs_;

    public:
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType  BlockMapperType;
      typedef SlaveDofs< DiscreteFunctionSpaceType, BlockMapperType > SlaveDofsType;


      // type of communication manager object which does communication
      typedef PetscDofMappings< SlaveDofsType >  PetscDofMappingType;
      typedef SingletonList< SlaveDofsType*, PetscDofMappingType > PetscDofMappingProviderType;

      explicit PetscSlaveDofProvider ( const DiscreteFunctionSpaceType &space )
      : BaseType( space ),
        dofMapping_( PetscDofMappingProviderType::getObject( slaveDofs_ ) )
      {
      }

      //! destructor
      ~PetscSlaveDofProvider()
      {
        // free mapping
        PetscDofMappingProviderType::removeObject( dofMapping_ );
      }

      const PetscDofMappingType& dofMapping() const { return dofMapping_; }
      PetscDofMappingType& dofMapping() { return dofMapping_; }

    protected:
      // the globally unique dof mapping
      PetscDofMappingType& dofMapping_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH
