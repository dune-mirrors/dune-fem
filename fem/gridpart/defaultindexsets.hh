#ifndef DUNE_FEM_DEFAULTINDEXSETS_HH
#define DUNE_FEM_DEFAULTINDEXSETS_HH

//- system includes 
#include <vector>
#include <rpc/rpc.h>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/alugrid/interfaces.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh> // for compatibility only
#include "../defaultindexsets.hh"

/** @file
 @author Robert Kloefkorn
 @brief Provides default index set implementations for Level- and
 LeafIndexsets used by ALUGrid. 
*/

namespace Dune {

//! compile time chooser for hierarchic or leaf index set
/** \deprecated */
template <class GridImp> 
class HierarchicIndexSetSelector
{
  
  // true if GridImp has HierarchicIndexSet 
  enum { hasHierarchicIndexSet = Conversion<GridImp,HasHierarchicIndexSet>::exists };

  template <class GridType, bool hasHSet> 
  struct HSetChooser
  {
    typedef typename GridType::Traits::LeafIndexSet IndexSetType;
    static const IndexSetType & hierarchicIndexSet(const GridType & grid)  
    { 
      return grid.leafIndexSet();
    }
  };
  
  template <class GridType> 
  struct HSetChooser<GridType,true>
  {
    typedef typename GridImp:: HierarchicIndexSet IndexSetType;
    static const IndexSetType & hierarchicIndexSet(const GridType & grid)  
    { 
      return grid.hierarchicIndexSet();
    }
  };
  
public: 
  //! \brief type of HierarchicIndexSet, default is LeafIndexSet 
  typedef typename HSetChooser<GridImp,hasHierarchicIndexSet>::IndexSetType HierarchicIndexSet; 
 
  //! \brief return reference to hierarchic index set 
  static const HierarchicIndexSet & hierarchicIndexSet(const GridImp & grid) 
  { 
    return HSetChooser<GridImp,hasHierarchicIndexSet>::hierarchicIndexSet(grid);
  }

  //! return true if index set can be used for adapitve calculations 
  static bool adaptive () { return hasHierarchicIndexSet; }
};

//! Wraps HierarchicIndex Sets of AlbertaGrid and ALUGrid 
/** \deprecated */
template <class GridType>
class WrappedHierarchicIndexSet
: public IndexSetWrapper< typename HierarchicIndexSetSelector<GridType> :: HierarchicIndexSet >
{
  // my type, to be revised 
  enum { myType = 0 };

  //! type of hset selector 
  typedef HierarchicIndexSetSelector<GridType> SelectorType;
  
  // my index set type 
  typedef typename SelectorType :: HierarchicIndexSet HSetType;

  typedef WrappedHierarchicIndexSet<GridType> ThisType;
public:
  //! number of codimensions 
  enum { ncodim = GridType::dimension + 1 };

  //! constructor 
  WrappedHierarchicIndexSet ( const GridType & grid , const int level =-1 ) 
    : IndexSetWrapper< HSetType > ( SelectorType :: hierarchicIndexSet(grid) 
                                  , SelectorType :: adaptive() ) {}
     
  //! return type (for Grape In/Output)
  static int type() { return myType; }

  //! returns reference to singleton 
  static ThisType & instance (const GridType & grid)
  { 
    static ThisType set(grid);
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
