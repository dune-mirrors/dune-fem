#ifndef DUNEFEM_INDEXSETS_HH
#define DUNEFEM_INDEXSETS_HH

//- system includes 
#include <iostream>
#include <string> 
#include <rpc/xdr.h>
#include <cassert>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/indexidset.hh>

//- Dune fem includes 
#include <dune/fem/gridpart/emptyindexset.hh>
#include <dune/fem/space/common/dofmanager.hh>

/** @file
 @brief Provides default index set class for persistent index sets. 
*/

namespace Dune {

/** \brief  class implementing the DUNE grid index set interface for any DUNE
            fem index set. However, one should not cast to this class as
            an interface class.
*/
template <class GridImp, class Imp, class IteratorTraits >
class DuneGridIndexSetAdapter 
  : public BartonNackmanInterface< 
        DuneGridIndexSetAdapter<GridImp, Imp, IteratorTraits > , 
        Imp > ,
    public EmptyIndexSet ,
    public IndexSet< GridImp, Imp, IteratorTraits > 
{
protected:
  typedef BartonNackmanInterface< 
            DuneGridIndexSetAdapter<GridImp, Imp, IteratorTraits > , 
            Imp > BaseType;
  using BaseType :: asImp;

  typedef IndexSet< GridImp, Imp, IteratorTraits > DuneIndexSetType;

  typedef DuneGridIndexSetAdapter<GridImp, Imp, IteratorTraits > ThisType;

public:
  //! type of index (i.e. unsigned int)
  typedef typename DuneIndexSetType :: IndexType IndexType;

  using DuneIndexSetType :: index; 

  friend class Conversion< ThisType, EmptyIndexSet> ;
  
  //! type of grid 
  typedef GridImp GridType;

  //! type of codimension 0 entity 
  typedef typename GridType :: template Codim<0> :: Entity  EntityCodim0Type; 

  //! constructor storing grid reference 
  DuneGridIndexSetAdapter(const GridType& grid) 
    : grid_(grid) 
  {}

  //! copy constructor 
  DuneGridIndexSetAdapter(const DuneGridIndexSetAdapter& other) 
    : grid_(other.grid_) 
  {}
  
protected:  
  //! reference to grid 
  const GridType& grid_;
  
public:
  //****************************************************************
  //
  //  INTERFACE METHODS for DUNE INDEX SETS 
  //
  //****************************************************************
  //! return global index of entity 
  template <class EntityType>
  inline IndexType index (const EntityType & en) const
  {
    // return index of entity 
    enum { codim = EntityType::codimension };
    return this->template index<codim> (en, 0);
  }

  //! return subIndex of given entity
  template <int cd>
  inline IndexType subIndex (const EntityCodim0Type & en, const int localNum) const
  {
    // return sub index of entity 
    return this->template index<cd> (en, localNum);
  }

  //////////////////////////////////////////////////////////////////
  //
  //  DUNE fem index method implementer interface  
  //
  //////////////////////////////////////////////////////////////////
protected:
  //! return index for entity  
  template <int codim, class EntityType>
  inline IndexType indexImp (const EntityType & entity , const int localNum) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().template indexImp<codim> ( entity , localNum ) );
    return asImp().template indexImp<codim> ( entity , localNum );
  } 

public:  
  //////////////////////////////////////////////////////////////////
  //
  //  DUNE fem index method user interface  
  //
  //////////////////////////////////////////////////////////////////
  //! return index for entity  
  template <int codim, class EntityType>
  inline IndexType index (const EntityType & entity , const int localNum) const
  {
    return this->template indexImp<codim> ( entity , localNum );
  } 

  //! insert new index for entity to set 
  inline void insertEntity(const EntityCodim0Type& entity )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().insertEntity( entity ) );
  }

  //! remove index for entity from index set 
  inline void removeEntity(const EntityCodim0Type& entity )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().removeEntity( entity ) );
  }
  
  // deprecated method  
  bool adaptive() const DUNE_DEPRECATED { return asImp().persistent(); } 
  // deprecated method  
  bool needsCompress () const DUNE_DEPRECATED { return asImp().consecutive(); }
};

/** \brief ConsecutivePersistentIndexSet is the base class for 
    all index sets that are persistent. Implementations of this type
    are for example all ConsecutivePersistenIndexSets. 
*/
template <class GridImp, class Imp, class IteratorTraits >
class PersistentIndexSet : 
  public DuneGridIndexSetAdapter<GridImp, Imp, IteratorTraits>
{
  // no copying 
  PersistentIndexSet(const PersistentIndexSet&);

  typedef Imp ImplementationType;
public:  
  //! type of base class 
  typedef DuneGridIndexSetAdapter<GridImp, Imp, IteratorTraits> BaseType;

  //! type of entity with codimension 0
  typedef typename BaseType :: EntityCodim0Type EntityCodim0Type;
  
  //! type of grid 
  typedef GridImp GridType;
  //! type of DoF manager
  typedef DofManager< GridType > DofManagerType;
  //! type of DoF manager factory
  typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

protected:
  using BaseType :: asImp ;

  // reference to dof manager 
  DofManagerType& dofManager_;

public:
  //! Conschdrugdor 
  inline PersistentIndexSet(const GridType & grid) 
    // here false, because methods have to be overloaded
    : BaseType(grid)
    , dofManager_( DofManagerFactoryType :: getDofManager( grid ) ) 
  {
    // add persistent index set to dofmanagers list 
    dofManager_.addIndexSet( asImp() );
  }

  //! destructor remoing index set from dof manager  
  inline ~PersistentIndexSet() 
  {
    // remove persistent index set from dofmanagers list 
    dofManager_.removeIndexSet( asImp() );
  }

  //! insert index for father, mark childs index for removal  
  inline void restrictLocal ( EntityCodim0Type &father, EntityCodim0Type &son, bool initialize )
  {
    // important, first remove old, because 
    // on father indices might be used aswell 
    asImp().removeEntity( son );
    asImp().insertEntity( father );
  }

  //! insert indices for children , mark fathers index for removal  
  inline void prolongLocal ( EntityCodim0Type &father, EntityCodim0Type &son, bool initialize )
  {
    // important, first remove old, because 
    // on children indices might be used aswell 
    asImp().removeEntity( father );
    asImp().insertEntity( son );
  }

  //! return true if the index set is persistent 
  inline bool persistent() const { return true; }
};


/** \brief ConsecutivePersistentIndexSet is the base class for all index sets that 
    are consecutive and also persistent. Implementations of this type
    are for example AdaptiveLeafIndexSet and DGAdaptiveLeafIndexSet. 
*/
template <class GridImp, class Imp, class IteratorTraits >
class ConsecutivePersistentIndexSet : public PersistentIndexSet<GridImp,Imp,IteratorTraits> 
{
  // no copying 
  ConsecutivePersistentIndexSet(const ConsecutivePersistentIndexSet&);

  typedef Imp ImplementationType;
public:  
  //! type of base class 
  typedef PersistentIndexSet<GridImp, Imp, IteratorTraits> BaseType;
  
  //! type of grid 
  typedef GridImp GridType;

  // use asImp from BaseType
  using BaseType :: asImp;
  
public:
  //! Conschdrugdor 
  inline ConsecutivePersistentIndexSet(const GridType & grid) 
    : BaseType(grid)
  {
  }

  //! returns true since we deal with a consecutive index set 
  inline bool consecutive () const { return true; }

  //! remove holes and make index set consecutive 
  inline bool compress() 
  {
    return asImp().compress();
  } 
};

} // end namespace Dune 
#endif
