#ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
#define DUNE_FEM_MANAGEDVECTORFUNCTION_HH

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
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

  protected:
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
      // allocate managed dof storage 
      std :: pair< DofStorageInterface *, DofVectorType * > memPair
        = allocateManagedDofStorage( dfSpace.grid(), dfSpace.mapper(), 
                                     name , (DofVectorType *) 0);
      memObject_ = memPair.first;
      return *(memPair.second);
    }
  };

}

#endif
