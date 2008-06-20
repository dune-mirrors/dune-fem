#ifndef DUNE_IDBASEDLEAFGRIDPART_HH
#define DUNE_IDBASEDLEAFGRIDPART_HH

//- Dune includes 
#include <dune/fem/gridpart/gridpart.hh>

//- local includes 
#include "idbasedleafindexset.hh"

namespace Dune {

/////////////////////////////////////////////////////////////////////////
//
//  --IdBasedLeafGridPart 
//
/////////////////////////////////////////////////////////////////////////

// forward deklaration of grid part 
template <class GridImp,PartitionIteratorType pitype>
struct IdBasedLeafGridPart;
template <class GridImp>
struct DefaultAdaptiveLeafIndexSet;

//! Type definitions for the LeafGridPart class
template <class GridImp,PartitionIteratorType pitype>
struct IdBasedLeafGridPartTraits {
  typedef GridImp GridType;
  typedef IdBasedLeafGridPart<GridImp,pitype> GridPartType;
  typedef IdBasedLeafIndexSet<GridImp> IndexSetType;

  typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;
  
  template <int cd>
  struct Codim {
    typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
  };

  //! \brief is true if grid on this view only has conformingi intersections 
  enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
};

/** @ingroup AdaptiveLeafGP
    \brief IdBasedLeafGridPart with
    indexset only for codimension 0 entities.

    Special implementation of the AdaptiveLeafGridPart with
    an underlying index set is only defined for 
    entities with codimension 0 for use with
    the Dune::DiscontinuousGalerkinSpace of FiniteVolumeSpaces.

    The underlying \ref IdBasedLeafIndexSet "index set" is
    a singleton for each different grid. 
    NOTE: The indices are stored in maps using the LocalIdSet of the grid 
          to generate these maps. 
    */
template <class GridImp, PartitionIteratorType pitype = Interior_Partition > 
class IdBasedLeafGridPart
: public GridPartDefault<IdBasedLeafGridPartTraits<GridImp,pitype> > 
{
public:
  //- Public typedefs and enums
  //! Type definitions
  typedef IdBasedLeafGridPartTraits<GridImp,pitype> Traits;
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

protected:
  // singleton provider 
  typedef SingletonList<const GridType* , IndexSetType > IndexSetProviderType;  

  // type of entity of codim 0
  typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

public:
  //- Public methods
  //! Constructor
  IdBasedLeafGridPart(GridType& grid) :
    GridPartDefault<Traits>(grid, IndexSetProviderType::getObject(&grid) )
  {}

  /** \brief Destrcutor removeing index set, if only one reference left, index set
      removed.  */
  ~IdBasedLeafGridPart() 
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
