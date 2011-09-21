#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH

//- system includes
#include <cassert>

//- dune-common includes 
#include <dune/common/bartonnackmanifcheck.hh>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>

//- dune-fem includes
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafindexset.hh>

namespace Dune
{
  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class, class, bool > struct FilteredGridPartTraits;

    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet = false >
    class FilteredGridPart;
   

    // FilteredGridPartIterator
    // ------------------------

    template< int codim, PartitionIteratorType pitype, class GridPartImp, class HostIteratorImp >
    class FilteredGridPartIterator
    {
      // type of this
      typedef FilteredGridPartIterator< codim, pitype, GridPartImp, HostIteratorImp > ThisType;

      // grid part type
      typedef GridPartImp GridPartType;

      // host iterator type
      typedef HostIteratorImp HostIteratorType;

      // entity pointer type
      typedef typename GridPartType::GridType::template Codim< codim >::EntityPointer EntityPointerType;

    public:
      // type of entity
      typedef typename HostIteratorType::Entity Entity;

      //! \brief constructor
      FilteredGridPartIterator ( const GridPartType & gridPart, const HostIteratorType & hostIterator )
      : gridPart_( gridPart ),
        hostIterator_( hostIterator ),
        hostEnd_( gridPart.hostGridPart().template end< codim, pitype >() )
      {
        if( done() ) 
          return;

        if( !gridPart.contains( *hostIterator_ ) )
          ++(*this);
      }

      //! \brief constructor
      FilteredGridPartIterator ( const ThisType & other )
      : gridPart_( other.gridPart_ ),
        hostIterator_( other.hostIterator_ ),
        hostEnd_( other.hostEnd_ )
      { }

      //! \brief assignment operator
      ThisType & operator= ( const ThisType & other )
      {
        assert( &gridPart_ == &other.gridPart_ );
        hostIterator_ = other.hostIterator_;
        hostEnd_ = other.hostEnd_;
        return *this;
      }

      //! \brief increment
      ThisType & operator++ ()
      {
        assert( !done() );
        do
        {
          ++hostIterator_;
        } while ( !done() && !contains() );
        return *this;
      }

      //! \brief return level
      int level () const
      {
        return hostIterator_.level();
      }

      const Entity & operator* () const
      {
        return *hostIterator_;
      }

      const Entity * operator-> () const
      {
        return &(*hostIterator_);
      }

      //! \brief cast to entity pointer
      operator EntityPointerType & ()
      {
        return hostIterator_;
      }

      //! \brief cast to const entity pointer
      operator const EntityPointerType & () const
      {
        return hostIterator_;
      }

      //! \brief check for equality
      bool operator== ( const ThisType & other ) const
      {
        return hostIterator_.operator==( other.hostIterator_ );
      }

      //! \brief check for inequality
      bool operator != ( const ThisType & other ) const
      {
        return !(*(this)==other);
      }

    private:
      // return true for end iterator
      bool done () const
      {
        return (hostIterator_ == hostEnd_ );
      }

      bool contains () const
      {
        assert( !done() );
        return gridPart().contains( *hostIterator_ );
      }

      // reference to grid part
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      const GridPartType & gridPart_;
      HostIteratorType hostIterator_;
      HostIteratorType hostEnd_;

    }; // end class FilteredGridPartIterator

    // IntersectionIteratorWrapper
    // ---------------------------

    template< class FilterType, class GridPartType, class HostIteratorType >
    class FilteredGridPartIntersectionIterator
    {
      // type of this
      typedef FilteredGridPartIntersectionIterator< FilterType, GridPartType, HostIteratorType > ThisType;

      // type of host intersecton
      typedef typename HostIteratorType::Intersection HostIntersection;

    public:
      //! \brief entity type
      typedef typename HostIntersection::Entity Entity;

      //! \brief entity pointer type
      typedef typename HostIntersection::EntityPointer EntityPointer;

      //! \brief geometry type
      typedef typename HostIntersection::Geometry Geometry;

      //! \brief local geometry type
      typedef typename HostIntersection::LocalGeometry LocalGeometry;

      //! \brief local coordinate type
      typedef typename HostIntersection::LocalCoordinate LocalCoordinate;

      //! \brief global coordinate type
      typedef typename HostIntersection::GlobalCoordinate GlobalCoordinate;

    protected:
      class NeighborInfo 
      {
        public:
        NeighborInfo ()
        : boundaryId_( -1 ),
          boundary_( false ),
          neighbor_(false) 
        { }
        
        NeighborInfo ( const NeighborInfo & org )
        : boundaryId_( org.boundaryId_ ),
          boundary_( org.boundary_ ),
          neighbor_( org.neighbor_ )
        { }
        
        NeighborInfo & operator = ( const NeighborInfo & org ) 
        {
          boundary_   = org.boundary_;
          boundaryId_ = org.boundaryId_; 
          neighbor_   = org.neighbor_;        
          return *this;
        }        

        int boundaryId_; 
        bool boundary_;
        bool neighbor_;        
      };

      // write information for current intersection 
      void writeNeighborInfo () 
      {
        if ( hostIterator()->neighbor() ) 
        { 
          if ( filter().interiorIntersection( *hostIterator() ) )
          {
            nInfo_.boundary_   = false;
            nInfo_.boundaryId_ = 0;
            nInfo_.neighbor_   = true;
          }
          else 
          {
            // otherwise get boundary information from filter 
            nInfo_.boundary_   = filter().intersectionBoundary( *hostIterator() );
            nInfo_.boundaryId_ = filter().intersectionBoundaryId( *hostIterator() );
            nInfo_.neighbor_   = filter().intersectionNeighbor( *hostIterator() );
          }
        }
        else 
        {
          // for real boundary get boundary from filter 
          nInfo_.boundary_   = true;
          nInfo_.boundaryId_ = filter().intersectionBoundaryId( *hostIterator() ); 
          nInfo_.neighbor_   = false;
        }    
      }


    public:
      //! \brief constructor 
      FilteredGridPartIntersectionIterator( const GridPartType & gridPart, 
                                            const HostIteratorType & hostIterator )
      : gridPart_( gridPart ),
        hostIterator_( hostIterator ),
        nInfo_()
      {
        if( !done() )
          writeNeighborInfo();
      }
        
      //! \brief copy constructor 
      FilteredGridPartIntersectionIterator( const ThisType & other )
      : gridPart_( other.gridPart_ ), 
        hostIterator_( other.hostIterator_ ),
        nInfo_( other.nInfo_ )
      { }
        
      //! \brief assignment operator 
      FilteredGridPartIntersectionIterator & operator = ( const ThisType & other ) 
      {
        gridPart_ = other.gridPart_;
        hostIterator_ = other.hostIterator_;
        nInfo_    = other.nInfo_; 
        return *this;
      }
        
      //! \brief increment intersection iterator 
      FilteredGridPartIntersectionIterator & operator++()
      { 
        assert( !done() );
        ++hostIterator_;
        if( !done() ) 
          writeNeighborInfo();
        return *this;
      }

      //! \brief check for equality 
      bool operator== ( const FilteredGridPartIntersectionIterator & other ) const
      {
        return hostIterator_.operator==( other.hostIterator_ );
      }

      //! \brief check for inequality 
      bool operator!= ( const FilteredGridPartIntersectionIterator & other ) const
      {
        return !(*this == other);
      }

      //! \brief overloaded boundary method 
      bool boundary () const
      {
        return nInfo_.boundary_;
      }

      //! \brief overloaded boundaryId method 
      int boundaryId () const
      {
        return nInfo_.boundaryId_;
      }

      //! \brief overloaded neighbor method 
      bool neighbor () const
      {
        return nInfo_.neighbor_;
      }

      //! \brief return inside entity
      EntityPointer inside () const
      {
        return hostIterator()->inside();
      }

      //! \brief return outside entity
      EntityPointer outside () const
      {
        return hostIterator()->outside();
      }

      //! \brief 
      bool conforming () const
      {
        return hostIterator()->conforming();
      }

      //! \brief return inside entity
      const LocalGeometry &geometryInInside () const
      {
        return hostIterator()->geometryInInside();
      }

      //! \brief return inside entity
      const LocalGeometry &geometryInOutside () const
      {
        return hostIterator()->geometryInOutside();
      }

      //! \brief return inside entity
      const Geometry &geometry () const
      {
        return hostIterator()->geometry();
      }

      //! \brief return inside entity
      GeometryType type () const
      {
        return hostIterator()->type();
      }

      //! \brief return inside entity
      int indexInInside () const
      {
        return hostIterator()->indexInInside();
      }

      //! \brief return inside entity
      int indexInOutside () const
      {
        return hostIterator()->indexInOutside();
      }

      //! \brief return inside entity
      GlobalCoordinate outerNormal ( const LocalCoordinate & local ) const
      {
        return hostIterator()->outerNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate & local ) const
      {
        return hostIterator()->integrationOuterNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate unitOuterNormal( const LocalCoordinate & local ) const
      {
        return hostIterator()->unitOuterNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate centerUnitOuterNormal( ) const
      {
        return hostIterator()->centerUnitOuterNormal( );
      }

      //! \brief type of Intersection 
      typedef ThisType Intersection;

      //! \brief dereference operator 
      const Intersection& operator *() const { return *this; }

      //! \brief de-pointer operator 
      const Intersection* operator ->() const { return this; }

    private:
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      bool done () const
      {
        return hostIterator_.operator == ( gridPart().hostGridPart().iend( *inside() ) );
      }

      // return reference to base class 
      HostIteratorType & hostIterator () 
      { 
        return hostIterator_;
      }

      // return reference to base class 
      const HostIteratorType & hostIterator () const
      {
        return hostIterator_;
      }

      const FilterType & filter () const
      {
        return gridPart().filter();
      }
      
      const GridPartType & gridPart_;
      HostIteratorType hostIterator_;
      NeighborInfo nInfo_;

    }; // end FilteredGridPartIntersectionIterator

    // FilteredGridPartIndexSetSelector
    // --------------------------------

    template < class FilteredGP, class HostGP, bool useFilteredIndexSet > 
    struct FilteredGridPartIndexSetSelector
    {
      typedef AdaptiveLeafIndexSet< FilteredGP > IndexSetType;

      static IndexSetType* create(const FilteredGP& gridPart) 
      {
        return new IndexSetType( gridPart );
      }

      template < class IndexSetPtr >
      static const IndexSetType & 
      indexSet ( const FilteredGP & gridPart, const IndexSetPtr * idxSetPtr )
      {
        assert( idxSetPtr );
        return *idxSetPtr;
      }
    };

    //! \brief when index set from gridpartimp is used return 0 
    template< class FilteredGP, class HostGP >
    struct FilteredGridPartIndexSetSelector< FilteredGP, HostGP, false >
    {
      typedef typename HostGP::IndexSetType IndexSetType;

      static IndexSetType* create(const FilteredGP& gridPart) 
      {
        return 0;
      }

      template < class IndexSetPtr >
      static const IndexSetType & 
      indexSet ( const FilteredGP & gridPart, const IndexSetPtr * )
      {
        return gridPart.hostGridPart().indexSet();
      }
    };


    // FilteredGridPartTraits
    // ----------------------

    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
    struct FilteredGridPartTraits
    {
      //! \brief type of grid part
      typedef FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > GridPartType;

      //! \brief grid part imp
      typedef HostGridPartImp HostGridPartType;

      //! \brief type of grid
      typedef typename HostGridPartType::GridType GridType;

      //! \brief export filter type
      typedef FilterImp FilterType;

      //! \brief type of entity
      typedef typename FilterType::EntityType EntityType;

      //! \brief index set use in this gridpart 
      typedef FilteredGridPartIndexSetSelector< GridPartType, HostGridPartType, useFilteredIndexSet > IndexSetSelectorType;

      //! \brief index set use in this gridpart 
      typedef typename IndexSetSelectorType::IndexSetType IndexSetType;
     
      //! \brief of host grid part intersection iterator type
      typedef typename HostGridPartType::Traits::IntersectionIteratorType HostIntersectionIteratorType;

      //! \brief type of intersection iterator 
      typedef FilteredGridPartIntersectionIterator< const FilterType, const GridPartType, HostIntersectionIteratorType > IntersectionIteratorType;

      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      //! \brief struct providing types of the iterators on codimension cd
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        class Partition
        {
          typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;

        public:
          typedef FilteredGridPartIterator< codim, pitype, const GridPartType, HostIteratorType > IteratorType;
        };

        typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;
      };

      //! \brief maximum partition type, the index set provides indices for
      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;

      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      //! \brief is true if grid on this view only has conforming intersections 
      static const bool conforming = HostGridPartType::conforming;
    };


    //***************************************************************************
    // 
    // FilteredGridPart
    //
    /** @addtogroup FilterGridPart 
     A FilteredGridPart is a subset of a GridPart and a GridPart itself. 
     The entities that belong to the FilteredGrid are defined by a 
     filter class. 
     
     Note that codim 0 entities have a method hasBoundaryIntersection().
     In general, this method will be inconsistent with the intersections
     returned by the filtered gridpart since entities are not wrapped. 
    **/

    /** @ingroup FilterGridPart
     @brief
     A FilteredGridPart allows to extract a set of entities from a grid
     satisfying a given constrainted defined through a filter class.
    **/ 


    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet > 
    class FilteredGridPart
    : public GridPartInterface< FilteredGridPartTraits< HostGridPartImp, FilterImp, useFilteredIndexSet > > 
    {
      // type of this
      typedef FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > ThisType;

    public:
      //- Public typedefs and enums    
      //! \brief traits class
      typedef FilteredGridPartTraits< HostGridPartImp, FilterImp, useFilteredIndexSet > Traits;
      
      //! \brief type of filter
      typedef FilterImp FilterType;

      // type of host grid part
      typedef typename Traits::HostGridPartType HostGridPartType;

      //! \brief grid type
      typedef typename Traits::GridType GridType;

      //! \brief index set type 
      typedef typename Traits::IndexSetType IndexSetType; 
      
      //! \brief intersection iterator type 
      typedef typename Traits:: IntersectionIteratorType IntersectionIteratorType;

      //! \brief intersection type
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      //! \brief grid view
      typedef GridView< GridPartViewTraits< ThisType > > GridViewType;

      //! \brief grid part typedefs 
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
        };
        typedef typename Partition< InteriorBorder_Partition > :: IteratorType IteratorType;
        typedef typename GridType::template Codim< codim >::Entity EntityType;
      };

    private:
      typedef typename Traits::IndexSetSelectorType IndexSetSelectorType;

      typedef typename Codim< 0 >::EntityType EntityType;

    public:
      //- Public methods
      //! \brief constructor
      FilteredGridPart ( HostGridPartType & hostGridPart, const FilterType & filter ) 
      : hostGridPart_( hostGridPart ),
        filter_( filter ),
        indexSetPtr_( 0 )
      {
        indexSetPtr_ = IndexSetSelectorType::create( *this );
      }

      //! \brief destructor 
      ~FilteredGridPart ()
      {
        if(  indexSetPtr_ )
          delete indexSetPtr_; 
      }

      //! \brief copy constructor
      FilteredGridPart ( const FilteredGridPart & other )
      : hostGridPart_( other.hostGridPart_ ), 
        filter_( other.filter_ ),
        indexSetPtr_ ( IndexSetSelectorType::create( *this ) )
      { }

      //! \brief return const reference to underlying grid
      const GridType & grid () const
      {
        return hostGridPart().grid();
      }

      //! \brief return reference to underlying grid
      GridType & grid ()
      {
        return hostGridPart().grid();
      }

      //! \brief return index set of this grid part 
      //         if IndexSetType is from host grid part the original index set is returned 
      const IndexSetType & indexSet() const 
      {
        return IndexSetSelectorType::indexSet( *this, indexSetPtr_ );
      } 
 
      //! \brief Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim >::IteratorType begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      //! \brief Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType begin () const
      {
        typedef typename Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
        return IteratorType( *this, hostGridPart().template begin< codim, pitype >() );
      }

      //! \brief Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim >::IteratorType end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! \brief End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType end () const
      {
        typedef typename Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
        return IteratorType( *this, hostGridPart().template end< codim, pitype >() );
      }

      //! \brief Returns maxlevel of the grid
      int level() const 
      { 
        return hostGridPart().level(); 
      }

      //! \brief ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType ibegin ( const EntityType & entity ) const 
      {
        return typename ThisType::IntersectionIteratorType( *this, hostGridPart().ibegin( entity ) );
      }
      
      //! \brief iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend ( const EntityType & entity ) const 
      {
        return typename ThisType::IntersectionIteratorType( *this, hostGridPart().iend( entity ) );
      }

      int boundaryId ( const IntersectionType & intersection ) const
      {
        return intersection.boundaryId();
      }

      //! \brief corresponding communication method for this grid part
      template < class DataHandleImp, class DataType >
      void communicate( CommDataHandleIF< DataHandleImp, DataType > & data, 
                        InterfaceType iftype, CommunicationDirection dir ) const 
      {
        this->grid().communicate( data, iftype, dir );
      }

      //! \brief return reference to filter 
      const FilterType & filter() const
      { 
        return filter_; 
      }

      //! \brief return reference to filter 
      FilterType & filter() 
      {
        return filter_; 
      }

      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        return filter().contains( entity );
      }
    
      HostGridPartType & hostGridPart ()
      {
        return hostGridPart_;
      }

      const HostGridPartType & hostGridPart () const
      {
        return hostGridPart_;
      }

    private: 
      HostGridPartType & hostGridPart_;
      const FilterType filter_;
      const IndexSetType * indexSetPtr_; 

    }; // end FilteredGridPart

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH
