#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTION_HH

#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartIntersection
    // ----------------------------

    template< class Filter, class HostIntersection >
    class FilteredGridPartIntersection
    {
      typedef FilteredGridPartIntersection< Filter, HostIntersection > ThisType;

    public:
      typedef Filter FilterType;
      typedef HostIntersection HostIntersectionType;

      static const int dimensionworld = HostIntersectionType::dimensionworld;
      static const int mydimension = HostIntersectionType::mydimension;

      typedef typename HostIntersectionType::ctype ctype;

      typedef typename HostIntersectionType::Entity Entity;
      typedef typename HostIntersectionType::Geometry Geometry;
      typedef typename HostIntersectionType::LocalGeometry LocalGeometry;

      typedef typename HostIntersectionType::LocalCoordinate LocalCoordinate;
      typedef typename HostIntersectionType::GlobalCoordinate GlobalCoordinate;

      FilteredGridPartIntersection () = default;

      FilteredGridPartIntersection ( const FilterType &filter, HostIntersectionType hostIntersection )
        : hostIntersection_( std::move( hostIntersection ) ),
          neighbor_( hostIntersection_.neighbor() ),
          boundary_( !neighbor_ )
      {
        if( neighbor_ )
        {
          if( !filter.interiorIntersection( hostIntersection_ ) )
          {
            neighbor_ = filter.intersectionNeighbor( hostIntersection_ );
            boundary_ = filter.intersectionBoundary( hostIntersection_ );
            boundaryId_ = filter.intersectionBoundaryId( hostIntersection_ );
          }
        }
        else
          boundaryId_ = filter.intersectionBoundaryId( hostIntersection_ );
      }

      bool equals ( const ThisType &other ) const { return (hostIntersection() == other.hostIntersection()); }

      bool boundary () const { return boundary_; }
      bool neighbor () const { return neighbor_; }

      int boundaryId () const { return boundaryId_; }

      std::size_t boundarySegmentIndex () const
      {
        DUNE_THROW( NotImplemented, "boundarySegmentIndex not implemented for FilteredGridPart, yet" );
      }

      Entity inside () const { return hostIntersection().inside(); }
      Entity outside () const { return hostIntersection().outside(); }

      bool conforming () const { return hostIntersection().conforming(); }

      LocalGeometry geometryInInside () const { return hostIntersection().geometryInInside(); }
      LocalGeometry geometryInOutside () const { return hostIntersection().geometryInOutside(); }

      Geometry geometry () const { return hostIntersection().geometry(); }
      GeometryType type () const { return hostIntersection().type(); }

      int indexInInside () const { return hostIntersection().indexInInside(); }
      int indexInOutside () const { return hostIntersection().indexInOutside(); }

      GlobalCoordinate outerNormal ( const LocalCoordinate & local ) const { return hostIntersection().outerNormal( local ); }
      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate & local ) const { return hostIntersection().integrationOuterNormal( local ); }
      GlobalCoordinate unitOuterNormal ( const LocalCoordinate & local ) const { return hostIntersection().unitOuterNormal( local ); }
      GlobalCoordinate centerUnitOuterNormal () const { return hostIntersection().centerUnitOuterNormal(); }

      const HostIntersectionType &hostIntersection () const { return hostIntersection_; }

    private:
      HostIntersectionType hostIntersection_;
      bool neighbor_ = false;
      bool boundary_ = false;
      int boundaryId_ = 0;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTION_HH
