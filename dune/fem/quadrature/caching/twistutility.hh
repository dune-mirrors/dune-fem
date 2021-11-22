#ifndef DUNE_FEM_TWISTUTILITY_HH
#define DUNE_FEM_TWISTUTILITY_HH

#include <cassert>

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/utility/hostgridaccess.hh>

// this also includes the forward declarations
#include <dune/fem/misc/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

   template< class Intersection >
    struct GridIntersectionAccess;

    template< class Grid, class IntersectionImpl >
    struct GridIntersectionAccess< Dune::Intersection< const Grid, IntersectionImpl > >
    {
      typedef Dune::Intersection< const Grid, IntersectionImpl > IntersectionType;
      typedef IntersectionType GridIntersectionType;

      static const GridIntersectionType &gridIntersection ( const IntersectionType &isec )
      {
        return isec;
      }
    };

    template< class Intersection >
    const typename GridIntersectionAccess< Intersection >::GridIntersectionType &
    gridIntersection ( const Intersection &isec )
    {
      return GridIntersectionAccess< Intersection >::gridIntersection( isec );
    }


    /** \brief TwistFreeTwistUtility provides the default implementation for twistfree grid
         such as Cartesian grids.
    */
    template< class Grid >
    struct TwistFreeTwistUtility
    {
      typedef Grid GridType;

      //! \brief return 0 for inner face
      template< class Intersection >
      static int twistInSelf ( const GridType &, const Intersection & )
      {
        return 0;
      }

      //! \brief return 0 for outer face
      template< class Intersection >
      static int twistInNeighbor ( const GridType &, const Intersection & )
      {
        return 0;
      }

      /** \brief return geometry type of inside or outside entity */
      template< class Intersection >
      static GeometryType elementGeometry ( const Intersection &intersection, const bool inside )
      {
        typedef Dune::Capabilities::hasSingleGeometryType< GridType > hasSingleGeometryType;
        if( hasSingleGeometryType::v )
          return GeometryType( hasSingleGeometryType::topologyId, GridType::dimension );
        else
          return inside ? intersection.inside().type() : intersection.outside().type();
      }
    };



    /** \brief Utility to get twist from IntersectionIterator,
        if provided by grid (i.e. AlbertaGrid, ALUGrid)
        otherwise return default values (correct for YASP).

        The twist (t) of a face is defined in the following way:
        - sign(t) gives information on the relationship between the
          orientation of the intersection geometry and the geometry
          of the corresponding codim 1 entity of the inside/outside
          entity:
          - sign(t)>=0: same orientation
          - sign(t)<0:  opposite orientation

        - The value of the twist gives information on the local numbering
          of the corners of the corresponding geometries. This value
          is only correctly defined for conforming grids, i.e.,
          the intersection is identical to an codim 1 entity of inside/outside.
          In this case we have the following definition:
          - sign(t)>=0: corner[0] of inside/outside face is equal to
                        corner[t] of intersection.
          - sign(t)<0:  corner[0] of inside/outside face is equal to
                        corner[t'] of intersection with t' = abs(t)+1.

        \note This class needs to be explicitly specialized for each grid!
    */
    template< class Grid >
    struct TwistUtility
     : public TwistFreeTwistUtility< Grid >
    {};

    // Specialization for AlbertaGrid
    // ------------------------------

#if HAVE_ALBERTA
    /** \brief Specialization of TwistUtility for AlbertaGrid.
    */
    template< int dim, int dimW >
    struct TwistUtility< AlbertaGrid< dim, dimW > >
    {
      typedef AlbertaGrid<dim, dimW> GridType;
      typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection  LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

      static const int dimension = GridType::dimension;

    public:
      //! \brief return twist for inner face
      static int twistInSelf ( const GridType &grid, const LeafIntersection &it )
      {
        return grid.getTwistInInside( it );
      }

      //! \brief return twist for outer face
      static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
      {
        return grid.getTwistInOutside( it );
      }

      /** \brief return element geometry type of inside or outside entity
      */
      template <class Intersection>
      static inline GeometryType
      elementGeometry(const Intersection& intersection,
                      const bool inside)
      {
        return Dune::GeometryTypes::simplex(dimension);
      }
    };
#endif // #if HAVE_ALBERTA



    // Specialization for ALUGrid
    // --------------------------

#if HAVE_DUNE_ALUGRID
    /** \brief Specialization of TwistUtility for ALUGrid.
    */
    template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct TwistUtility< ALUGrid< dim, dimw, elType, refineType, Comm > >
    {
      typedef ALUGrid< dim, dimw, elType, refineType, Comm > GridType;

    public:
      //! \brief return twist for inner face
      template< class Intersection >
      static inline int twistInSelf(const GridType & grid, const Intersection& isec)
      {
        const auto& intersection = gridIntersection( isec );
        assert( dim == 2 ? (intersection.impl().twistInInside() == 0 ||
                            intersection.impl().twistInInside() == 1 ) : true );
        return intersection.impl().twistInInside();
      }

      //! \brief return twist for outer face
      template< class Intersection >
      static inline int twistInNeighbor(const GridType &grid, const Intersection& isec )
      {
        const auto& intersection = gridIntersection( isec );
        assert( dim == 2 ? (intersection.impl().twistInOutside() == 0 ||
                            intersection.impl().twistInOutside() == 1 ) : true );
        return intersection.impl().twistInOutside();
      }

      /** \brief return element geometry type of inside or outside entity
      */
      template <class Intersection>
      static inline GeometryType
      elementGeometry(const Intersection& intersection,
                      const bool inside)
      {
        return GeometryType( Dune::Capabilities::hasSingleGeometryType< GridType > :: topologyId,
                             dim );
      }

    private:
      TwistUtility(const TwistUtility&);
      TwistUtility& operator=(const TwistUtility&);
    };
#endif // #if HAVE_DUNE_ALUGRID



    // Specialization for UGGrid
    // -------------------------

#if HAVE_UG
    template< int dim >
    struct TwistUtility< UGGrid< dim > >
    {
      typedef UGGrid< dim > GridType;

      typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

      static int twistInSelf ( const GridType &grid, const LeafIntersection &it );
      static int twistInSelf ( const GridType &grid, const LevelIntersection &it );

      static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it );
      static int twistInNeighbor ( const GridType &grid, const LevelIntersection &it );

      template< class Intersection >
      static GeometryType
      elementGeometry ( const Intersection &intersection, const bool inside )
      {
        return (inside ? intersection.inside().type() : intersection.outside().type());
      }
    };
#endif // #ifdef ENABLE_UG



    // Specialization for GeometryGrid
    // -------------------------------

    template< class HostGrid, class CoordFunction, class Allocator >
    struct TwistUtility< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Allocator > GridType;
      typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    private:
      typedef TwistUtility< HostGrid > HostTwistUtility;
      typedef Dune::HostGridAccess< GridType > HostGridAccess;

    public:
      //! \brief return twist for inner face
      template< class Intersection >
      static int twistInSelf ( const GridType &grid, const Intersection &intersection )
      {
        return HostTwistUtility::twistInSelf( grid.hostGrid(), HostGridAccess::hostIntersection( intersection ) );
      }

      //! \brief return twist for outer face
      template< class Intersection >
      static int twistInNeighbor ( const GridType &grid, const Intersection &intersection )
      {
        return HostTwistUtility::twistInNeighbor( grid.hostGrid(), HostGridAccess::hostIntersection( intersection ) );
      }

      /** \brief return element geometry type of inside or outside entity */
      template< class Intersection >
      static GeometryType elementGeometry ( const Intersection &intersection, bool inside )
      {
        return HostTwistUtility::elementGeometry( HostGridAccess::hostIntersection( intersection ), inside );
      }
    };

  }  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TWISTUTILITY_HH
