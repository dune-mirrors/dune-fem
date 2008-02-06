#ifndef DUNE_PERSISTENTINDEXSETS_HH
#define DUNE_PERSISTENTINDEXSETS_HH

//- Dune includes 
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/fem/space/common/dofmanager.hh>

/** @file
 @author Robert Kloefkorn
 @brief Provides default index set class for persistent index sets. 
*/

namespace Dune {

/*!
  The DefaultGridIndexSet is a wrapper for the grid index which can be 
  index of entity or globalIndex. The DofMapper uses an IndexSet for
  mapping the dofs, so we can hide the real grid index behind the index
  set. Furthermore if a grid doesn't provide the consecutive index set
  then this can be calculated in the IndexSet. These two following index
  sets are just the identiy to the grid indices. 

  The DefaultGridIndexSetBase defines some methods that are needed for
  index sets that cope with adaptation, but aren't needed for the following 
  index set, so most of this methods do notin'.
*/

//! default base class for index set implementations for FR numerics
template <class GridImp, class Imp>
class PersistentIndexSet : public DefaultGridIndexSetBase<GridImp>
{
  // no copying 
  PersistentIndexSet(const PersistentIndexSet&);

  typedef Imp ImplementationType;
public:  
  //! type of base class 
  typedef DefaultGridIndexSetBase<GridImp> BaseType;
  
  //! type of grid 
  typedef GridImp GridType;
  //! type of DoF manager
  typedef DofManager< GridType > DofManagerType;
  //! type of DoF manager factory
  typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

protected:
  // reference to dof manager 
  DofManagerType& dofManager_;

  //! cast to implementation class for add to dofmanager 
  //! otherwise dofmanager will use the wrong methods from the 
  //! default empty base class 
  ImplementationType& asImp() 
  { 
    return static_cast<ImplementationType &> (*this); 
  }
  
public:
  //! Conschdrugdor 
  PersistentIndexSet(const GridType & grid) 
    // here false, because methods have to be overloaded
    : BaseType(grid)
    , dofManager_( DofManagerFactoryType :: getDofManager( grid ) ) 
  {
    // add persistent index set to dofmanagers list 
    dofManager_.addIndexSet( asImp() );
  }

  ~PersistentIndexSet() 
  {
    // remove persistent index set from dofmanagers list 
    dofManager_.removeIndexSet( asImp() );
  }
};

} // end namespace Dune 
#endif
