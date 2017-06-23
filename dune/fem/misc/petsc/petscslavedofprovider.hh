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
    class PetscSlaveDofProvider
    {
    public:
      typedef DFSpace DiscreteFunctionSpaceType;

      typedef PetscSlaveDofProvider< DiscreteFunctionSpaceType > ThisType;

      typedef typename DiscreteFunctionSpaceType :: BlockMapperType  BlockMapperType;

      // type of communication manager object which does communication
      typedef PetscDofMappings<DiscreteFunctionSpaceType> PetscDofMappingType;
      typedef SingletonList< const DiscreteFunctionSpaceType*, PetscDofMappingType > PetscDofMappingProviderType;

      explicit PetscSlaveDofProvider ( const DiscreteFunctionSpaceType &space )
      : dofMapping_( PetscDofMappingProviderType::getObject( &space ) )
      {
        dofMapping_.update( space );
      }

      //! destructor
      ~PetscSlaveDofProvider()
      {
        // free mapping
        PetscDofMappingProviderType::removeObject( dofMapping_ );
      }

      //! update dof mapping
      void update(const DiscreteFunctionSpaceType &space)
      {
        dofMapping_.update(space);
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
