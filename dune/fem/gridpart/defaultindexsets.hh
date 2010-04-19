#ifndef DUNE_FEM_DEFAULTINDEXSETS_HH
#define DUNE_FEM_DEFAULTINDEXSETS_HH

//- system includes 
#include <vector>
#include <rpc/rpc.h>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh> // for compatibility only
#include <dune/grid/alugrid/defaultindexsets.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/mpimanager.hh>

/** @file
 @author Robert Kloefkorn
 @brief Provides default index set implementations for Level- and
 LeafIndexsets used by ALUGrid. 
*/

namespace Dune {

//! Wraps the interface methods of indexsets and adds the addiotnal needed
//! functions 
template <class IndexSetImp> 
class IndexSetWrapper : public DefaultEmptyIndexSet 
{
public:
  //! The types of the iterator 
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename IndexSetImp::template Codim<cd>::
        template Partition<pitype>::Iterator  Iterator;
    };
  };

  //! store const reference to set 
  IndexSetWrapper(const IndexSetImp & set, bool adaptive = false) 
    : DefaultEmptyIndexSet(adaptive)
    , set_(set) 
  {}
  
  //! store const reference to set 
  IndexSetWrapper(const IndexSetWrapper<IndexSetImp> & s) 
    : DefaultEmptyIndexSet(s.adaptive_), set_(s.set_) {}

  //! return persistent status 
  bool persistent () const { return false; }
 
  //! return size of set for codim  
  int size ( GeometryType type ) const   
  {
    return set_.size(type);
  }

  //! return size of grid entities per level and codim 
  int size ( int codim ) const   
  {
    return set_.size(codim);
  }

  //! return index of en 
  template <class EntityType> 
  int index (const EntityType & en) const
  {
    return set_.index(en); 
  }

  //! return sub index of given entities sub entity with codim and number 
  template< int codim, class EntityType >
  int DUNE_DEPRECATED subIndex ( const EntityType &entity, int num ) const
  {
    return set_.template subIndex< codim >( entity, num );
  }

  template< class Entity >
  int subIndex ( const Entity &entity, int num, unsigned int codim ) const
  {
    return set_.subIndex( entity, num, codim );
  }

  //! wrap geomTypes method of set 
  const std::vector< GeometryType > & geomTypes (int codim) const 
  {
    return set_.geomTypes(codim); 
  }

  //! returns true if this set provides an index for given entity
  template<class EntityType>
  bool contains (const EntityType& en) const
  {
    return set_.contains(en); 
  }

  //! return index 
  template <int codim, class EntityType> 
  int index (const EntityType & en, int num) const
  {
    enum { enCodim = EntityType::codimension };
    return IndexWrapper<IndexSetImp,EntityType,enCodim,codim>::index(set_,en,num);
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin () const
  { 
    return set_.template begin<cd,pitype> ();
  }

  /** @brief Iterator to one past the last entity of given codim for partition type
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end () const
  { 
    return set_.template end<cd,pitype> ();
  }

private: 
  const IndexSetImp & set_; 
};
//! Wraps LevelIndexSet for use with LagrangeFunctionSpace 
template <class GridType>
class WrappedLevelIndexSet 
: public IndexSetWrapper< typename GridType :: Traits:: LevelIndexSet > 
{
  // my type, to be revised 
  enum { myType = 1 };

  typedef typename GridType :: Traits :: LevelIndexSet LevelIndexSetType;
  typedef WrappedLevelIndexSet < GridType >  ThisType;
public:
  
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };

  //! Constructor getting grid and level for Index Set 
  WrappedLevelIndexSet (const GridType & grid , const int level ) 
    : IndexSetWrapper<  LevelIndexSetType > (grid.levelIndexSet(level)) {}
 
  //! return type of index set (for input/output)
  static int type() { return myType; }
  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid,grid.maxLevel());
    return set;
  }

};

//! compile time chooser for hierarchic or leaf index set
/** \deprecated */
template< class Grid >
class HierarchicIndexSetSelector
{
  template< bool >
  struct HierarchicIndexSetGetter
  {
    typedef typename Grid::HierarchicIndexSet IndexSet;

    static const IndexSet &indexSet ( const Grid &grid )
    {
      return grid.hierarchicIndexSet();
    }
  };

  template< bool >
  struct LeafIndexSetGetter
  {
    typedef typename Grid::LeafIndexSet IndexSet;

    static const IndexSet &indexSet ( const Grid &grid )
    {
      if( Dune::MPIManager::rank() == 0 )
        std::cerr << "Warning: Grid does not provide a HierarchicIndexSet, using LeafIndexSet instead." << std::endl;
      return grid.leafIndexSet();
    }
  };

  static const bool hasHierarchicIndexSet = Capabilities::hasHierarchicIndexSet< Grid >::v;
  typedef typename SelectType< hasHierarchicIndexSet, HierarchicIndexSetGetter< true >, LeafIndexSetGetter< false > >::Type IndexSetGetter;
  
public: 
  //! \brief type of HierarchicIndexSet, default is LeafIndexSet
  typedef typename IndexSetGetter::IndexSet HierarchicIndexSet;
 
  //! \brief return reference to hierarchic index set 
  static const HierarchicIndexSet &hierarchicIndexSet ( const Grid &grid )
  { 
    return IndexSetGetter::indexSet( grid );
  }

  //! return true if index set can be used for adapitve calculations 
  static bool adaptive ()
  {
    return hasHierarchicIndexSet;
  }
};



//! Wraps HierarchicIndex Sets of AlbertaGrid and ALUGrid 
/** \deprecated */
template< class GridType >
class WrappedHierarchicIndexSet
: public IndexSetWrapper< typename HierarchicIndexSetSelector< GridType >::HierarchicIndexSet >
{
  typedef WrappedHierarchicIndexSet< GridType > ThisType;
  typedef IndexSetWrapper< typename HierarchicIndexSetSelector< GridType >::HierarchicIndexSet > BaseType;

  // my type, to be revised 
  enum { myType = 0 };

  // type of hierarchic index set selector
  typedef HierarchicIndexSetSelector< GridType > SelectorType;
  
public:
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };

  //! constructor 
  WrappedHierarchicIndexSet ( const GridType &grid, const int level =-1 )
  : BaseType( SelectorType::hierarchicIndexSet( grid ), SelectorType::adaptive() )
  {}
     
  //! return type (for Grape In/Output)
  static int type ()
  {
    return myType;
  }

  //! returns reference to singleton
  static ThisType &instance ( const GridType &grid )
  {
    static ThisType set( grid );
    std::cerr << "Warning: WrappedHierarchicIndexSet::instance( grid ) can only handle one grid." << std::endl;
    return set;
  }
};



//! Wraps LeafIndexSet of Dune Grids for use with LagrangeFunctionSpace 
template <class GridType>
class WrappedLeafIndexSet
:  public IndexSetWrapper<typename GridType :: Traits :: LeafIndexSet> 
{
  // my type, to be revised 
  enum { myType = 5 };

  // my index set type 
  typedef typename GridType :: Traits :: LeafIndexSet IndexSetType;
  typedef WrappedLeafIndexSet<GridType> ThisType;
public:
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };
  //! constructor 
  WrappedLeafIndexSet ( const GridType & grid , const int level =-1 ) 
    : IndexSetWrapper < IndexSetType > (grid.leafIndexSet()) {}

  //! constructor taking grid part 
  template <class GridPartType> 
  WrappedLeafIndexSet ( const GridPartType& gridPart ) 
    : IndexSetWrapper < IndexSetType > (gridPart.grid().leafIndexSet()) {}

  //! return type (for Grape In/Output)
  static int type() { return myType; }
  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid,grid.maxLevel());
    return set;
  }
};

} // end namespace Dune 

#endif
