#ifndef DUNE_FEM_PERIODICGRIDPART_HH
#define DUNE_FEM_PERIODICGRIDPART_HH

#include <dune/grid/common/grid.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/periodicindexset.hh>

namespace Dune
{
  
  template< class Grid >
  class PeriodicLeafGridPart;



  template< class Grid >
  class PeriodicLeafIntersectionIterator
  {
    typedef PeriodicLeafIntersectionIterator< Grid > ThisType;
      
    friend class PeriodicLeafGridPart< Grid >;

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
    explicit
    PeriodicLeafIntersectionIterator ( WrappedIteratorType wrappedIterator )
    : wrappedIterator_( wrappedIterator )
    {}

  public:
    PeriodicLeafIntersectionIterator ( const ThisType &other )
    : wrappedIterator_( other.wrappedIterator_ )
    {}

    ThisType &operator= ( const ThisType &other )
    {
      wrappedIterator_ = other.wrappedIterator_;
      return *this;
    }

    ThisType &operator++ ()
    {
      ++wrappedIterator_;
      return *this;
    }

    const Intersection &operator* () const
    {
      return *this;
    }

    const Intersection *operator-> () const
    {
      return this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return (wrappedIterator_ == other.wrappedIterator_);
    }

    bool operator!= ( const ThisType &other ) const
    {
      return (wrappedIterator_ != other.wrappedIterator_);
    }

    bool boundary () const
    {
      // IntersectionIterator specifies that this should return true at the boundary
      return false;
    }

    int boundaryId () const
    {
      return 0;
    }

    int neighbor () const
    {
      // IntersectionIterator specifies that we should return true!
      return wrappedIterator_->neighbor();
    }

    EntityPointer inside () const
    {
      return wrappedIterator_->inside();
    }

    EntityPointer outside () const
    {
      if( wrappedIterator_->neighbor() )
        return wrappedIterator_->outside();
      else
        DUNE_THROW( NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                    "outside on boundary not implemented yet." );
    }

    const LocalGeometry &intersectionSelfLocal () const
    {
      return wrappedIterator_->intersectionSelfLocal();
    }

    const LocalGeometry &intersectionNeighborLocal () const
    {
      if( wrappedIterator_->neighbor() )
        return wrappedIterator_->intersectionNeighborLocal();
      else
        DUNE_THROW( NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                    "outside on boundary not implemented yet." );
    }

    const Geometry &intersectionGlobal () const
    {
      return wrappedIterator_->intersectionGlobal();
    }

    int numberInSelf () const
    {
      return wrappedIterator_->numberInSelf();
    }

    int numberInNeighbor () const
    {
      return wrappedIterator_->numberInNeighbor();
    }

    DomainType outerNormal ( const LocalDomainType &x ) const
    {
      return wrappedIterator_->outerNormal( x );
    }

    DomainType integrationOuterNormal ( const LocalDomainType &x ) const
    {
      return wrappedIterator_->integrationOuterNormal( x );
    }

    DomainType unitOuterNormal ( const LocalDomainType &x ) const
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

  
          
  template< class Grid >
  class PeriodicLeafGridPartTraits
  {
    typedef PeriodicLeafGridPartTraits< Grid > ThisType;

  public:
    typedef Grid GridType;

    typedef PeriodicLeafGridPart< GridType > GridPartType;
   
    typedef PeriodicLeafIndexSet< GridType > IndexSetType;

    static const PartitionIteratorType indexSetPartitionType = All_Partition;

    typedef PeriodicLeafIntersectionIterator< GridType >
      IntersectionIteratorType;

    template< int codim >
    struct Codim
    {
      // type of the entity (should be wrapped, too)
      typedef typename GridType :: template Codim< codim > :: Entity EntityType;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename GridType
          :: template Codim< codim > :: template Partition< pitype > :: LeafIterator
          IteratorType;
      };
    };

    //! is the grid partition conforming?
    static const bool conforming = Capabilities :: isLeafwiseConforming< GridType > :: v;
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
  template< class Grid >
  class PeriodicLeafGridPart
  : public GridPartDefault< PeriodicLeafGridPartTraits< Grid > >
  {
    typedef PeriodicLeafGridPart< Grid > ThisType;
    typedef GridPartDefault< PeriodicLeafGridPartTraits< Grid > > BaseType;

  public:
    //! type of traits
    typedef typename BaseType :: Traits Traits;
    
    //! type of the underlying grid
    typedef typename Traits :: GridType GridType;

  public:
    //! type of the index set
    typedef typename Traits :: IndexSetType IndexSetType;

    //! type of intersection iterators
    typedef typename Traits :: IntersectionIteratorType IntersectionIteratorType;
    
    template< int codim >
    struct Codim
    : public BaseType :: template Codim< codim >
    {
      typedef typename Traits :: template Codim< codim > :: EntityType EntityType;
    };

  protected:
    IndexSetType indexSet_;

  public:
    //! Constructor wrapping a grid in the grid partition
    PeriodicLeafGridPart ( GridType &grid )
    : BaseType( grid, indexSet_ ),
      indexSet_( grid )
    {}

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    begin () const
    {
      return BaseType :: template begin< codim >();
    }

    //! Begin iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    begin () const
    {
      return (*this).grid().template leafbegin< codim, pitype >();
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    end () const
    {
      return BaseType :: template end< codim >();
    }

    //! End iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    end () const
    {
      return (*this).grid().template leafend< codim, pitype >();
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
    IntersectionIteratorType
    ibegin ( const typename Codim< 0 > :: EntityType &entity ) const
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
    IntersectionIteratorType
    iend ( const typename Codim< 0 > :: EntityType &entity ) const
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
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      this->grid().communicate( data, iftype, dir );
    }
  };

}

#endif
