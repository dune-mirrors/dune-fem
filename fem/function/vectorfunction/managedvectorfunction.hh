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
    typedef DofManager< GridType > DofManagerType;

  protected:
    DofManagerType *dofManager_;
    MemObjectInterface *memObject_;

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
      assert( memObject_ );
      dofManager_->removeDofSet( *memObject_ );
    }

    inline void enableDofCompression ()
    {
      assert( memObject_ );
      memObject_->enableDofCompression();
    }

  private:
    inline DofVectorType &
    allocDofVector ( const std :: string &name,
                     const DiscreteFunctionSpaceType &dfSpace )
    {
      dofManager_ = &DofManagerFactory< DofManagerType >
                       :: getDofManager( dfSpace.grid() );

      std :: pair< MemObjectInterface *, DofVectorType * > memPair
        = dofManager_->addDofSet( (DofVector *)0, dfSpace.mapper(), name );
      memObject_ = memPair.first;
      return *(memPair.second);
    }
  };

}

#endif
