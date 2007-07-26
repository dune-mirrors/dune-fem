#ifndef DUNE_FEM_PERIODICGRIDPART_HH
#define DUNE_FEM_PERIODICGRIDPART_HH

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridpart.hh>

#include "periodicindexset.hh"

namespace Dune
{
  
  template< class GridImp, PartitionIteratorType pitype >
  class PeriodicLeafGridPart;

  
          
  template< class GridImp, PartitionIteratorType pitype >
  class PeriodicLeafGridPartTraits
  {
  public:
    //! type of the underlying grid
    typedef GridImp GridType;

  private:
    typedef PeriodicLeafGridPartTraits< GridType, pitype > ThisType;

  public:
    //! type of the grid partition (Barton-Nackman)
    typedef PeriodicLeafGridPart< GridType, pitype > GridPartType;
   
    //! type of the index set
    typedef WrappedPeriodicLeafIndexSet< GridType > IndexSetType;

    //! type of codim 0 entities
    typedef typename GridType :: template Codim< 0 > :: Entity Codim0EntityType;

    //! type of intersection iterators
    typedef typename Codim0EntityType :: LeafIntersectionIterator IntersectionIteratorType;

    template< int codim >
    struct Codim
    {
      typedef typename GridType
        :: template Codim< codim > :: template Partition< pitype > :: LeafIterator
        IteratorType;
    };

    //! is the grid partition conforming?
    enum { conforming = Capabilities :: isLeafwiseConforming< GridType > :: v };
  };
          


  /*! \class PeriodicGridPart
   *  \brief A grid partition identifying opposite faces of the unit cube
   *
   *  Using any grid of the unit cube, this grid partition makes the grid
   *  periodic by identifying the indices of subentities of opposite faces.
   *
   *  \note Since the underlying grid does not know about the periodic boundary,
   *        refinement may break conformity (global refinement should work,
   *        though).
   *
   *  \todo Return correct neighbors for entities with boundary intersections.
   */
  template< class GridImp, PartitionIteratorType pitype = Interior_Partition >
  class PeriodicLeafGridPart
  : public GridPartDefault< PeriodicLeafGridPartTraits< GridImp, pitype > >
  {
  public:
    //! type of traits
    typedef PeriodicLeafGridPartTraits< GridImp, pitype > TraitsType;
    
    //! type of the underlying grid
    typedef typename TraitsType :: GridType GridType;

  private:
    typedef PeriodicLeafGridPart< GridType, pitype > ThisType;
    typedef GridPartDefault< TraitsType > BaseType;

  public:
    //! type of the index set
    typedef typename TraitsType :: IndexSetType IndexSetType;

    //! type of codim 0 entities
    typedef typename TraitsType :: Codim0EntityType Codim0EntityType;

    //! type of intersection iterators
    typedef typename TraitsType :: IntersectionIteratorType IntersectionIteratorType;
    
    template< int codim >
    struct Codim
    {
      typedef typename TraitsType :: template Codim< codim > :: IteratorType IteratorType;
    };

  protected:
    IndexSetType indexSet_;

  public:
    //! Constructor wrapping a grid in the grid partition
    PeriodicLeafGridPart ( GridType &grid )
    : BaseType( grid, indexSet_ ),
      indexSet_( grid )
    {
    }

    //! Begin iterator for periodic grid
    template< int codim >
    typename Codim< codim > :: IteratorType begin () const
    {
      return this->grid().template leafbegin< codim, pitype >();
    }

    //! End iterator for periodic grid
    template< int codim >
    typename Codim< codim > :: IteratorType end () const
    {
      return this->grid().template leafend< codim, pitype >();
    }

    //! Begin intersection iterator for an entity
    IntersectionIteratorType ibegin ( const Codim0EntityType &entity ) const
    {
      return entity.ileafbegin();
    }

    //! End intersection iterator for an entity
    IntersectionIteratorType iend ( const Codim0EntityType &entity ) const
    {
      return entity.ileafend();
    }

    //! Deliver maximum level of grid
    int level () const
    {
      return this->grid().maxLevel();
    }

    //! Communication Method for this grid partition
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data, InterfaceType iftype, CommunicationDirection dir ) const
    {
      this->grid().communicate( data, iftype, dir );
    }
  };

}

#endif
