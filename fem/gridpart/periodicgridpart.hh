#ifndef DUNE_FEM_PERIODICGRIDPART_HH
#define DUNE_FEM_PERIODICGRIDPART_HH

#include <dune/grid/common/grid.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/periodicindexset.hh>

namespace Dune
{
  
  template< class GridImp, PartitionIteratorType pitype >
  class PeriodicLeafGridPart;



  template< class Grid >
  class PeriodicLeafIntersectionIterator
  {
    typedef PeriodicLeafIntersectionIterator< Grid > ThisType;
      
    template< class, PartitionIteratorType >
    friend class PeriodicLeafGridPart;

  public:
    typedef Grid GridType;

    typedef ThisType Intersection;

  protected:
    typedef typename GridType :: template Codim< 0 > :: Entity Codim0EntityType;
   
    typedef typename Codim0EntityType :: LeafIntersectionIterator
      WrappedIteratorType;

  public:
    typedef typename WrappedIteratorType :: ctype ctype;

    enum
    {
      dimension = WrappedIteratorType :: dimension,
      dimensionworld = WrappedIteratorType :: dimensionworld
    };

    typedef FieldVector< ctype, dimensionworld > DomainType;

    typedef FieldVector< ctype, dimensionworld - 1 > LocalDomainType;
    
    typedef typename WrappedIteratorType :: Entity Entity;

    typedef typename WrappedIteratorType :: EntityPointer EntityPointer;

    typedef typename WrappedIteratorType :: Geometry Geometry;

    typedef typename WrappedIteratorType :: LocalGeometry LocalGeometry;

    typedef typename WrappedIteratorType :: ImplementationType
      ImplementationType;

  protected:
    WrappedIteratorType wrappedIterator_;

  protected:
    inline explicit
    PeriodicLeafIntersectionIterator ( WrappedIteratorType wrappedIterator )
    : wrappedIterator_( wrappedIterator )
    {}

  public:
    inline PeriodicLeafIntersectionIterator ( const ThisType &other )
    : wrappedIterator_( other.wrappedIterator_ )
    {}

    inline ThisType &operator= ( const ThisType &other )
    {
      wrappedIterator_ = other.wrappedIterator_;
      return *this;
    }

    inline ThisType &operator++ ()
    {
      ++wrappedIterator_;
      return *this;
    }

    inline const Intersection &operator* () const
    {
      return *this;
    }

    inline const Intersection *operator-> () const
    {
      return this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return (wrappedIterator_ == other.wrappedIterator_);
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return (wrappedIterator_ != other.wrappedIterator_);
    }

    inline bool boundary () const
    {
      // IntersectionIterator specifies that this should return true at the boundary
      return false;
    }

    inline int boundaryId () const
    {
      return 0;
    }

    inline int neighbor () const
    {
      // IntersectionIterator specifies that we should return true!
      return wrappedIterator_->neighbor();
    }

    inline EntityPointer inside () const
    {
      return wrappedIterator_->inside();
    }

    inline EntityPointer outside () const
    {
      if( wrappedIterator_->neighbor() )
        return wrappedIterator_->outside();
      else
        DUNE_THROW( NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                    "outside on boundary not implemented yet." );
    }

    inline const LocalGeometry &intersectionSelfLocal () const
    {
      return wrappedIterator_->intersectionSelfLocal();
    }

    inline const LocalGeometry &intersectionNeighborLocal () const
    {
      if( wrappedIterator_->neighbor() )
        return wrappedIterator_->intersectionNeighborLocal();
      else
        DUNE_THROW( NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                    "outside on boundary not implemented yet." );
    }

    inline const Geometry &intersectionGlobal () const
    {
      return wrappedIterator_->intersectionGlobal();
    }

    inline int numberInSelf () const
    {
      return wrappedIterator_->numberInSelf();
    }

    inline int numberInNeighbor () const
    {
      return wrappedIterator_->numberInNeighbor();
    }

    inline DomainType outerNormal ( const LocalDomainType &x ) const
    {
      return wrappedIterator_->outerNormal( x );
    }

    inline DomainType integrationOuterNormal ( const LocalDomainType &x ) const
    {
      return wrappedIterator_->integrationOuterNormal( x );
    }

    inline DomainType unitOuterNormal ( const LocalDomainType &x ) const
    {
      return wrappedIterator_->unitOuterNormal( x );
    }

  protected:
    const ImplementationType &getRealImp () const
    {
      return wrappedIterator_.getRealImp();
    }

    ImplementationType &getRealImp ()
    {
      return wrappedIterator_.getRealImp();
    }
  };

  
          
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
    //typedef typename Codim0EntityType :: LeafIntersectionIterator
    //  IntersectionIteratorType;
    typedef PeriodicLeafIntersectionIterator< GridType >
      IntersectionIteratorType;

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
          


  /*! \addtogroup PeriodicGridPart
   *  \class PeriodicLeafGridPart
   *  \brief A grid partition identifying opposite faces of the unit cube
   *
   *  Using any grid of the unit cube, this grid partition makes the grid
   *  periodic by identifying the indices of subentities of opposite faces.
   *
   *  \note Since the underlying grid does not know about the periodic boundary,
   *        refinement may break conformity (global refinement should work,
   *        though).
   *
   *  \note This grid partition says that there is no boundary. In DUNE, however,
   *        periodic boundaries shall be implemented as boundaries with ghost
   *        entities (however the FEM codes usually only check if an intersection
   *        is a boundary.
   *
   *  \todo The entity needs also to be wrapped, so that hasBoundaryIntersections
   *        always returns false
   *
   *  \todo Return correct neighbors for entities with boundary intersections.
   *
   *  \newimplementation Allows to construct globally refined grids for the
   *                     unitcube with periodic boundaries.
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

    //! type of codim 0 entities (these must be wrapped, too)
    typedef typename TraitsType :: Codim0EntityType Codim0EntityType;

    //! type of intersection iterators
    typedef typename TraitsType :: IntersectionIteratorType IntersectionIteratorType;
    
    template< int codim >
    struct Codim
    {
      typedef typename TraitsType :: template Codim< codim > :: IteratorType
        IteratorType;
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

    /** \brief begin intersection iterator for an entity
     *
     *  \note The intersection iterators always return boundary = false
     *        and neighbor = true.
     * 
     *  \param[in]  entity  entity the intersection iterator is requested for
     *
     *  \returns a begin intersection iterator
     */
    IntersectionIteratorType ibegin ( const Codim0EntityType &entity ) const
    {
      return IntersectionIteratorType( entity.ileafbegin() );
    }
    
    /** \brief end intersection iterator for an entity
     *
     *  \note The intersection iterators always return boundary = false
     *        and neighbor = true.
     * 
     *  \param[in]  entity  entity the intersection iterator is requested for
     *
     *  \returns an end intersection iterator
     */
    IntersectionIteratorType iend ( const Codim0EntityType &entity ) const
    {
      return IntersectionIteratorType( entity.ileafend() );
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
