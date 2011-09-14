#ifndef DUNE_FILTEREDGRID_HH
#define DUNE_FILTEREDGRID_HH

//- system includes
#include <vector>
#include <cassert>

//- dune-common includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/typetraits.hh>

// dune-grid includes
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

// dune-fem includes
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafindexset.hh>

namespace Dune
{
  namespace Fem
  {

    // forward declarations
    // --------------------

    template < class, class > struct DefaultFilterTraits;
    template< class > class FilterInterface;
    template< class > class FilterDefaultImplementation;
    template< class, bool > struct FilteredGridPartTraits;

    // filtered grid part
    // ------------------

    template< class FilterImp, bool useFilteredIndexSet = false >
    class FilteredGridPart;


    // DefaultFilterTraits
    // -------------------

    //! \brief type definitions 
    template < class FilterImp, class GridPartImp >
    struct DefaultFilterTraits 
    {
      //! \brief filter type
      typedef FilterImp FilterType;

      //! \brief grid part type
      typedef GridPartImp GridPartType; 

      template < int cd >
      struct Codim 
      {
        //! \brief entity type
        typedef typename GridPartType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief entity type for codim 0 
      typedef typename Codim< 0 >::EntityType EntityType;   

    }; // end DefaultFilterTraits



    // FilterInterface
    // ---------------

    /** \ingroup FilterGridPart
     *  \brief   Interface class for filter to use with a Dune::FilteredGridPart
     */
    template< class FilterTraits >
    class FilterInterface
    {
      typedef FilterInterface< FilterTraits > ThisType;

      friend class FilterDefaultImplementation< FilterTraits >;

    public:
      //! \brief Type of the filter implementation
      typedef typename FilterTraits :: FilterType FilterType;

      //! \brief type of original grid part 
      typedef typename FilterTraits :: GridPartType GridPartType;

      template< int cd >
      struct Codim
      {
        typedef typename FilterTraits::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of entity with codim=0
      typedef typename Codim< 0 >::EntityType EntityType;

    private: 
      // constructor
      FilterInterface () { }

      // copy constructor
      FilterInterface ( const ThisType & );

      // assignment operator
      ThisType & operator= ( const ThisType & );

    public:
      //! \brief returns true if the given entity of the pointer in the domain 
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().contains( entity ) );
        return asImp().contains( entity );
      }

      //! \brief returns true if the given entity of the pointer in the domain 
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        enum { cc = Entity::codimension };
        CHECK_INTERFACE_IMPLEMENTATION( asImp().contains< cc >( entity ) );
        return asImp().contains< cc >( entity );
      }
      
      //! \brief returns true if an intersection is interior 
      //         (allows boundarys within a given domain)
      template< class Intersection >
      bool interiorIntersection ( const Intersection &intersection ) const
      {
        return asImp().interiorIntersection( intersection );
      }

      //! \brief returns true if an intersection is a boundary intersection 
      template< class Intersection >
      bool intersectionBoundary( const Intersection &intersection ) const
      {
        return asImp().intersectionBoundary( intersection );
      }
      
      //! \brief returns the boundary id for an intersection 
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection &intersection ) const
      {
        return asImp().intersectionBoundaryId( intersection );
      }

      //! \brief returns true if for an intersection a neighbor exsits 
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection &intersection ) const
      {
        return asImp().intersectionNeighbor( intersection );
      }

      //! \brief return object instance of filter for object creation  
      static FilterType createObject( const GridPartType &gridPart )
      {
        return FilterType :: createObject( gridPart );
      }

    protected:
      FilterType &asImp ()
      {
        return static_cast< FilterType & >( *this );
      }
      
      const FilterType &asImp () const
      {
        return static_cast< const FilterType & >( *this );
      }
    };


    // FilterDefaultImplementation
    // ---------------------------

    template< class FilterTraits >
    class FilterDefaultImplementation
    : public FilterInterface< FilterTraits >
    {
      typedef FilterDefaultImplementation< FilterTraits > ThisType;
      typedef FilterInterface< FilterTraits > BaseType;

    public:
      //! \brief type of the filter implementation
      typedef typename BaseType::FilterType FilterType;

      //! \brief type of original grid part 
      typedef typename BaseType::GridPartType GridPartType;
     
      template< int cd >
      struct Codim
      {
        typedef typename BaseType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of codim 0 entity 
      typedef typename BaseType::EntityType EntityType;
       
    protected:
      // constructor
      FilterDefaultImplementation () { }

      using BaseType :: asImp;

      // copy constructor
      FilterDefaultImplementation ( const ThisType & );

      // assignment operator 
      ThisType &operator= ( const ThisType & );

    public:
      using BaseType::contains;

      //! \brief default implementation returns contains from neighbor
      template< class Intersection >
      bool interiorIntersection( const Intersection &intersection ) const
      {
        typedef typename GridPartType::GridType::template Codim< 0 >::EntityPointer EntityPointerType;
        const EntityPointerType outside = intersection.outside();
        return asImp().contains( *outside );
      }

      //! \brief default createObject method calling FilterType( gridPart )
      static FilterType createObject( const GridPartType & gridPart )
      {
        return FilterType( gridPart );
      }
     
      //! \brief returns true if the given entity of the pointer in the domain 
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & ) const;

      //! \brief returns true if an intersection is a boundary intersection 
      template< class Intersection >
      bool intersectionBoundary( const Intersection & ) const;
     
      //! \brief returns the boundary id for an intersection 
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection & ) const;

      //! \brief returns true if for an intersection a neighbor exsits 
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & ) const;
    };

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
        hostIterator_( hostIterator )
      {
        if( done() ) 
          return;

        if( !gridPart.contains( *hostIterator_ ) )
          ++(*this);
      }

      //! \brief increment
      ThisType & operator++ ()
      {
        assert( !done() );
        ++hostIterator_;
        if( done() )
          return *this;
        if( !gridPart().contains( *hostIterator_ ) )
          ++(*this);
        return *this;
      }

      //! \brief return level
      int level () const
      {
        return hostIterator_.level();
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
        return hostIterator_.operator==( gridPart().hostGridPart().template end< codim, pitype >() );
      }

      // reference to grid part
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      const GridPartType & gridPart_;
      HostIteratorType hostIterator_;

    }; // end class FilteredGridPartIterator

    // IntersectionIteratorWrapper
    // ---------------------------

    template< class FilterType, class GridPartType, class HostIteratorType >
    class FilteredGridPartIntersectionIterator
    : public HostIteratorType
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
        inline NeighborInfo ()
        : boundaryId_( -1 ),
          boundary_( false ),
          neighbor_(false) 
        { }
        
        inline NeighborInfo ( const NeighborInfo & org )
        : boundaryId_( org.boundaryId_ ),
          boundary_( org.boundary_ ),
          neighbor_( org.neighbor_ )
        { }
        
        inline NeighborInfo & operator = ( const NeighborInfo & org ) 
        {
          boundary_   = org.boundary_;
          boundaryId_ = org.boundaryId_; 
          neighbor_   = org.neighbor_;        
          return *this;
        }        

        int boundaryId_; 
        bool boundary_;
        bool neighbor_;        
      } nInfo_;

      // write information for current intersection 
      inline void writeNeighborInfo() 
      {
        if ( asBase()->neighbor() ) 
        { 
          if ( filter().interiorIntersection( *asBase() ) )
          {
            nInfo_.boundary_   = false;
            nInfo_.boundaryId_ = 0;
            nInfo_.neighbor_   = true;
          }
          else 
          {
            // otherwise get boundary information from filter 
            nInfo_.boundary_   = filter().intersectionBoundary( *asBase() );
            nInfo_.boundaryId_ = filter().intersectionBoundaryId( *asBase() );
            nInfo_.neighbor_   = filter().intersectionNeighbor( *asBase() );
          }
        }
        else 
        {
          // for real boundary get boundary from filter 
          nInfo_.boundary_   = true;
          nInfo_.boundaryId_ = filter().intersectionBoundaryId( *asBase() );
          nInfo_.neighbor_   = false;
        }    
      }


    public:
      //! \brief constructor 
      FilteredGridPartIntersectionIterator( const GridPartType & gridPart, 
                                            const HostIteratorType & hostIterator )
      : HostIteratorType( hostIterator ),
        nInfo_(),
        gridPart_( gridPart )
      {
        if( !done() )
          writeNeighborInfo();
      }
        
      //! \brief copy constructor 
      FilteredGridPartIntersectionIterator( const ThisType & other )
      : HostIteratorType( other ),
        nInfo_( other.nInfo_ ),
        gridPart_( other.gridPart_ )
      { }
        
      //! \brief assignment operator 
      FilteredGridPartIntersectionIterator & operator = ( const ThisType & other ) 
      {
        HostIteratorType::operator = ( other );
        nInfo_    = other.nInfo_; 
        gridPart_ = other.gridPart_;
        return *this;
      }
        
      //! \brief increment intersection iterator 
      inline FilteredGridPartIntersectionIterator & operator++()
      { 
        assert( !done() );
        asBase().operator++();
        if( done() ) 
          return *this; 
        writeNeighborInfo();
        return *this;
      }

      //! \brief check for equality 
      bool operator== ( const FilteredGridPartIntersectionIterator & other ) const
      {
        return asBase().operator==( other.asBase() );
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
        return asBase()->inside();
      }

      //! \brief return outside entity
      EntityPointer outside () const
      {
        return asBase()->outside();
      }

      //! \brief 
      bool conforming () const
      {
        return asBase()->conforming();
      }

      //! \brief return inside entity
      const LocalGeometry &geometryInInside () const
      {
        return asBase()->geometryInInside();
      }

      //! \brief return inside entity
      const LocalGeometry &geometryInOutside () const
      {
        return asBase()->geometryInOutside();
      }

      //! \brief return inside entity
      const Geometry &geometry () const
      {
        return asBase()->geometry();
      }

      //! \brief return inside entity
      GeometryType type () const
      {
        return asBase()->type();
      }

      //! \brief return inside entity
      int indexInInside () const
      {
        return asBase()->indexInInside();
      }

      //! \brief return inside entity
      int indexInOutside () const
      {
        return asBase()->indexInOutside();
      }

      //! \brief return inside entity
      GlobalCoordinate outerNormal ( const LocalCoordinate & local ) const
      {
        return asBase()->outerNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate & local ) const
      {
        return asBase()->integrationOuterNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate unitOuterNormal( const LocalCoordinate & local ) const
      {
        return asBase()->unitOuterNormal( local );
      }

      //! \brief return inside entity
      GlobalCoordinate centerUnitOuterNormal( ) const
      {
        return asBase()->centerUnitOuterNormal( );
      }

      //! \brief type of Intersection 
      typedef ThisType Intersection;

      //! \brief dereference operator 
      inline const Intersection& operator *() const { return *this; }

      //! \brief de-pointer operator 
      inline const Intersection* operator ->() const { return this; }

    private:
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      bool done () const
      {
        return asBase().operator == ( gridPart().hostGridPart().iend( *inside() ) );
      }

      // return reference to base class 
      inline HostIteratorType & asBase() 
      { 
        return static_cast< HostIteratorType & >( *this ); 
      }

      // return reference to base class 
      inline const HostIteratorType & asBase () const
      {
        return static_cast< const HostIteratorType & >( *this );
      }

      const FilterType & filter () const
      {
        return gridPart().filter();
      }
      
      const GridPartType & gridPart_;

    }; // end FilteredGridPartIntersectionIterator

    // FilteredGridPartIndexSetSelector
    // --------------------------------

    template < class HostGridPart, bool useFilteredIndexSet > 
    struct FilteredGridPartIndexSetSelector
    {
      typedef AdaptiveLeafIndexSet< HostGridPart > IndexSetType;

      static IndexSetType* create (const HostGridPart & gridPart) 
      {
        return new IndexSetType( gridPart );
      }

      template < class IndexSetPtr >
      static inline const IndexSetType & 
      indexSet ( const HostGridPart & gridPart, const IndexSetPtr * idxSetPtr )
      {
        assert( idxSetPtr );
        return *idxSetPtr;
      }
    };

    //! \brief when index set from gridpartimp is used return 0 
    template< class HostGridPart >
    struct FilteredGridPartIndexSetSelector< HostGridPart, false >
    {
      typedef typename HostGridPart::IndexSetType IndexSetType;

      static IndexSetType* create ( const HostGridPart & ) 
      {
        return 0;
      }

      template < class IndexSetPtr >
      static inline const IndexSetType & 
      indexSet ( const HostGridPart & gridPart, const IndexSetPtr * )
      {
        return gridPart.indexSet();
      }
    };


    // FilteredGridPartTraits
    // ----------------------

    template< class FilterImp, bool useFilteredIndexSet >
    struct FilteredGridPartTraits
    {
      //! \brief type of grid part
      typedef FilteredGridPart< FilterImp, useFilteredIndexSet > GridPartType;

      //! \brief grid part imp
      typedef typename FilterImp::GridPartType HostGridPartType;

      //! \brief type of grid
      typedef typename HostGridPartType::GridType GridType;

      //! \brief export filter type
      typedef FilterImp FilterType;

      //! \brief type of entity
      typedef typename FilterType::EntityType EntityType;

      //! \brief index set use in this gridpart 
      typedef FilteredGridPartIndexSetSelector< HostGridPartType, useFilteredIndexSet > IndexSetSelectorType;

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
     A FilteredGridPart is a subset of a GridPart and a GridPart itself 
     (without iterators for codim \f$\neq\f$ 0). 
     The codim 0 entities that belong to the FilteredGrid are defined by a 
     filter class. On a codim 0 entitiy there is a method 
       hasBoundaryIntersection().
     This method will not work correctly since the entity is not wrapped. 
    **/

    /** @ingroup FilterGridPart
     @brief
     A FilteredGridPart allows to extract a set of entities from a grid
     satisfying a given constrainted defined through a filter class.
    **/ 


    template< class FilterImp, bool useFilteredIndexSet > 
    class FilteredGridPart
    : public GridPartInterface< FilteredGridPartTraits< FilterImp, useFilteredIndexSet > > 
    {
      // type of this
      typedef FilteredGridPart< FilterImp, useFilteredIndexSet > ThisType;

    public:
      //- Public typedefs and enums    
      //! \brief traits class
      typedef FilteredGridPartTraits< FilterImp, useFilteredIndexSet > Traits;
      
      //! \brief type of filter
      typedef FilterImp FilterType;

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

      // type of host grid part
      typedef typename Traits::HostGridPartType HostGridPartType;

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
        indexSetPtr_ = IndexSetSelectorType::create( hostGridPart_ );
      }

      //! \brief destructor 
      ~FilteredGridPart ()
      {
         if( indexSetPtr_ )
          delete indexSetPtr_; 
      }

      //! \brief copy constructor
      FilteredGridPart ( const FilteredGridPart & other )
      : hostGridPart_( other.hostGridPart_ ), 
        filter_( other.filter_ ),
        indexSetPtr_( IndexSetSelectorType::create( hostGridPart_ ) )
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
        return IndexSetSelectorType::indexSet( hostGridPart(), indexSetPtr_ );
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
      inline IntersectionIteratorType ibegin ( const EntityType & entity ) const 
      {
        return typename ThisType::IntersectionIteratorType( *this, hostGridPart().ibegin( entity ) );
      }
      
      //! \brief iend of corresponding intersection iterator for given entity
      inline IntersectionIteratorType iend ( const EntityType & entity ) const 
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
      const FilterType & filter_;
      const IndexSetType * indexSetPtr_; 

    }; // end FilteredGridPart

  }  // end namespace Fem

}  // end namespace Dune

#endif
