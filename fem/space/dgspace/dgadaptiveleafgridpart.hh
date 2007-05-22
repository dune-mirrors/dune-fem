#ifndef DGADAPTIVELEAFGRIDPART_HH
#define DGADAPTIVELEAFGRIDPART_HH

//- Dune includes 
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/idbasedleafindexset.hh>

//- local includes 
#include "dgadaptiveleafindexset.hh"

namespace Dune { 

// forward deklaration 
template <class GridImp,PartitionIteratorType pitype>
struct DGAdaptiveLeafGridPart;

//! default is to use DGAdaptiveLeafIndexSet
template <class GridImp, bool unStructured, bool isGood> 
struct GoodGridChooser
{
  // choose the adative index based on hierarhic index set
  typedef DGAdaptiveLeafIndexSet<GridImp> IndexSetType;
};

// the same for shitty unstructured grids 
template <class GridImp> 
struct GoodGridChooser<GridImp,true,false>
{
  // choose the leaf index set based on local ids 
  typedef IdBasedLeafIndexSet<GridImp> IndexSetType;
};

// for structured grids with no hierarchic index set choose 
// the grids leaf index set 
template <class GridImp> 
struct GoodGridChooser<GridImp,false,false>
{
#if HAVE_MPI
  // choose the adative index based on leaf index set
  typedef DGAdaptiveLeafIndexSet<GridImp> IndexSetType;
#else 
  // the grids leaf index set wrapper for good 
  typedef WrappedLeafIndexSet<GridImp> IndexSetType;
#endif
};

template <class GridType>
class DGAdaptiveLeafIndexSet; 

//! Type definitions for the LeafGridPart class
template <class GridImp,PartitionIteratorType pitype>
struct DGAdaptiveLeafGridPartTraits {
  
  typedef GridImp GridType;
  
  typedef DGAdaptiveLeafGridPart<GridImp,pitype> GridPartType;
  // choose index set dependend on grid type  
  typedef GoodGridChooser<GridImp,
           Capabilities::IsUnstructured<GridType>::v,
           Conversion<GridImp,HasHierarchicIndexSet>::exists
         > IndexSetChooserType;
              
  typedef typename IndexSetChooserType :: IndexSetType IndexSetType;

  typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;
  
  template <int cd>
  struct Codim {
    typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
  };
      
  //! \brief is true if grid on this view only has conforming intersections 
  enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
};
/** @ingroup AdaptiveLeafGP
    @brief
    Special implementation of the AdaptiveLeafGridPart with
    an underlying index set only defined for entities with codimension 0.
    The underlying index set is
    a singleton for each different grid. */
template <class GridImp, PartitionIteratorType pitype = Interior_Partition > 
class DGAdaptiveLeafGridPart
: public GridPartDefault<DGAdaptiveLeafGridPartTraits<GridImp,pitype> > 
{
  //! singleton list , key type is const pointer to grid 
  typedef SingletonList< const GridImp* ,
          typename DGAdaptiveLeafGridPartTraits<GridImp,pitype>::IndexSetType > IndexSetProviderType;  
public:
  //- Public typedefs and enums
  //! Type definitions
  typedef DGAdaptiveLeafGridPartTraits<GridImp,pitype> Traits;
  //! Grid implementation type
  typedef typename Traits::GridType GridType;
  //! The leaf index set of the grid implementation
  typedef typename Traits::IndexSetType IndexSetType;
  
  //! The corresponding IntersectionIterator 
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;
  
  //! Struct providing types of the leaf iterators on codimension cd
  template <int cd>
  struct Codim {
    typedef typename Traits::template Codim<cd>::IteratorType IteratorType;
  };

  //! \brief is true if grid on this view only has conforming intersections 
  enum { conforming = Traits :: conforming };

private:
  typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

public:
  //- Public methods
  //! Constructor
  DGAdaptiveLeafGridPart(GridType& grid) :
    GridPartDefault<Traits>(grid, IndexSetProviderType::getObject(&grid) )
  {}

  /** \brief Desctrutor removing index set. When only one reference is
      left, index set object is deleted. */
  ~DGAdaptiveLeafGridPart() 
  { 
    IndexSetProviderType::removeObject(this->indexSet());
  }

  //! Begin iterator on the leaf level
  template <int cd>
  typename Traits::template Codim<cd>::IteratorType begin() const {
    return this->grid().template leafbegin<cd,pitype>();
  }

  //! End iterator on the leaf level
  template <int cd>
  typename Traits::template Codim<cd>::IteratorType end() const {
    return this->grid().template leafend<cd,pitype>();
  }

  //! ibegin of corresponding intersection iterator for given entity
  IntersectionIteratorType ibegin(const EntityCodim0Type & en) const 
  {
    return en.ileafbegin();
  }
  
  //! iend of corresponding intersection iterator for given entity
  IntersectionIteratorType iend(const EntityCodim0Type & en) const 
  {
    return en.ileafend();
  }

  //! Returns maxlevel of the grid
  int level() const { return this->grid().maxLevel(); }

  //! corresponding communication method for this grid part
  template <class DataHandleImp,class DataType>
  void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                   InterfaceType iftype, CommunicationDirection dir) const 
  {
    this->grid().communicate(data,iftype,dir);
  }
};

}// end namespace Dune
#endif
