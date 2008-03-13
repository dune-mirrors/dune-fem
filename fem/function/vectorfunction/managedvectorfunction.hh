#ifndef DUNE_FEM_MANAGEDVECTORFUNCTION_HH
#define DUNE_FEM_MANAGEDVECTORFUNCTION_HH

#include <dune/fem/function/vectorfunction/vectorfunction.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  template< class DiscreteFunction >
  class ManagedDiscreteFunction;


  template< class DiscreteFunctionSpace, class DofVector >
  class ManagedDiscreteFunction
    < VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
  : public VectorDiscreteFunction< DiscreteFunctionSpace, DofVector >
  {
    typedef VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > BaseType;
    typedef ManagedDiscreteFunction< BaseType > ThisType;

  public:
    typedef typename BaseType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename BaseType :: DofVectorType DofVectorType;

  protected:
    typedef typename DiscreteFunctionSpaceType :: GridPartType :: GridType
      GridType;
    typedef DofManager< GridType > DofManagerType;

  protected:
    DofManagerType *dofManager_;
    std :: pair< MemObjectInterface *, DofVectorType * > memPair_;

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
      if( memPair_.first )
        dofManager_->removeDofSet( *(memPair_.first) );
    }

  private:
    inline DofVectorType &
    allocDofVector ( const std :: string &name,
                     const DiscreteFunctionSpaceType &dfSpace )
    {
      dofManager_ = &DofManagerFactory< DofManagerType >
                       :: getDofManager( dfSpace.grid() );
      memPair_
        = dofManager_->addDofSet( (DofVector *)0, dfSpace.mapper(), name );
      return *(memPair_.second);
    }
  };

}

#endif
