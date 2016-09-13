#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH

//- system includes
#include <cassert>

//- dune-grid includes
#include <dune/grid/common/intersectioniterator.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartIntersectionIterator
    // ------------------------------------

    template< class FilterType, class GridPartType, class HostIteratorType >
    class FilteredGridPartIntersectionIterator
    {
      // type of this
      typedef FilteredGridPartIntersectionIterator< FilterType, GridPartType, HostIteratorType > ThisType;

      // type of host intersecton
      typedef typename HostIteratorType::Intersection HostIntersection;

    public:
      //! \brief dimension
      static const int dimension = HostIntersection::dimension;
      //! \brief world dimension
      static const int dimensionworld = HostIntersection::dimensionworld;
      static const int mydimension = dimension - 1;

      //! \brief single coordinate type
      typedef typename HostIntersection::ctype ctype;

      //! \brief entity type
      typedef typename HostIntersection::Entity Entity;

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
                                            const Entity &en,
                                            const HostIteratorType & hostIterator )
      : gridPart_( gridPart ),
        hostIterator_( hostIterator ),
        endIterator_( gridPart.hostGridPart().iend( en ) ),
        nInfo_()
      {
        if( !done() )
          writeNeighborInfo();
      }

      //! \brief copy constructor
      FilteredGridPartIntersectionIterator( const ThisType & other )
      : gridPart_( other.gridPart_ ),
        hostIterator_( other.hostIterator_ ),
        endIterator_( other.endIterator_ ),
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

      std::size_t boundarySegmentIndex () const
      {
        DUNE_THROW( NotImplemented, "boundarySegmentIndex not implemented for FilteredGridPart, yet" );
      }

      //! \brief overloaded neighbor method
      bool neighbor () const
      {
        return nInfo_.neighbor_;
      }

      //! \brief return inside entity
      Entity inside () const
      {
        return hostIterator()->inside();
      }

      //! \brief return outside entity
      Entity outside () const
      {
        return hostIterator()->outside();
      }

      //! \brief
      bool conforming () const
      {
        return hostIterator()->conforming();
      }

      //! \brief return inside entity
      LocalGeometry geometryInInside () const
      {
        return hostIterator()->geometryInInside();
      }

      //! \brief return inside entity
      LocalGeometry geometryInOutside () const
      {
        return hostIterator()->geometryInOutside();
      }

      //! \brief return inside entity
      Geometry geometry () const
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

      // impl method is needed for MetaTwistUtility
      const Intersection& impl() const
      {
        return *this;
      }

      // hostIntersection method is needed for MetaTwistUtility
      HostIntersection hostIntersection() const
      {
        return *hostIterator();
      }

    private:
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      bool done () const
      {
        return hostIterator_.operator ==( endIterator_ );
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
      HostIteratorType hostIterator_, endIterator_;
      NeighborInfo nInfo_;

    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH
