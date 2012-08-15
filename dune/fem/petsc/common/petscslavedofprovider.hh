// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH
#define DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH

#include <dune/fem/function/common/scalarproducts.hh>


#if defined HAVE_PETSC

#include <dune/fem/petsc/common/petsccommon.hh>
#include <dune/fem/petsc/common/petscdofmappings.hh>


namespace Dune {
namespace Fem {

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
  typedef Dune::SlaveDofs< DiscreteFunctionSpaceType, BlockMapperType > SlaveDofsType;


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

#endif // #if defined HAVE_PETSC

#endif // DUNE_FEM_PETSCSLAVEDOFPROVIDER_HH