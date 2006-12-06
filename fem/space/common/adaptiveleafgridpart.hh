#ifndef DUNE_ADAPTIVELEAFGRIDPART_HH
#define DUNE_ADAPTIVELEAFGRIDPART_HH

//- Dune includes 
#include <dune/grid/common/gridpart.hh>

//- local includes 
#include "adaptiveleafindexset.hh"

namespace Dune {

/////////////////////////////////////////////////////////////////////////
//
//  --AdaptiveLeafIndexGridPart 
//
/////////////////////////////////////////////////////////////////////////

// forward deklaration of grid part 
template <class GridImp,PartitionIteratorType pitype>
struct AdaptiveLeafGridPart;
template <class GridImp>
struct AdaptiveLeafIndexSet;

//! Type definitions for the LeafGridPart class
template <class GridImp,PartitionIteratorType pitype>
struct AdaptiveLeafGridPartTraits {
  typedef GridImp GridType;
  typedef AdaptiveLeafGridPart<GridImp,pitype> GridPartType;
  typedef AdaptiveLeafIndexSet<GridImp> IndexSetType;

  typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;
  
  template <int cd>
  struct Codim {
    typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
  };

  //! \brief is true if grid on this view only has conforming intersections 
  enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
};

/** \brief GridPart for AdaptiveLeafIndexSet. Used underlying index set is
    singleton for each grid object.   
*/
template <class GridImp, PartitionIteratorType pitype = Interior_Partition > 
class AdaptiveLeafGridPart
: public GridPartDefault<AdaptiveLeafGridPartTraits<GridImp,pitype> > 
{
  // singleton list for index set , key type is const pointer to grid 
  typedef SingletonList<const GridImp *,
          typename AdaptiveLeafGridPartTraits<GridImp,pitype>::IndexSetType > 
            IndexSetProviderType;  
public:
  //- Public typedefs and enums
  //! Type definitions
  typedef AdaptiveLeafGridPartTraits<GridImp,pitype> Traits;
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
  AdaptiveLeafGridPart(GridType& grid) :
    GridPartDefault<Traits>(grid, IndexSetProviderType::getObject(&grid) )
  {}

  /** \brief Destrcutor removeing index set, if only one reference left, index set
      removed.  */
  ~AdaptiveLeafGridPart() 
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

} // end namespace Dune 
#endif
