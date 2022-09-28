#ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
#define DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/declaration.hh>
#endif // #if HAVE_DUNE_SPGRID

#if HAVE_OPM_GRID
#include <opm/grid/polyhedralgrid/declaration.hh>
#endif // #if HAVE_OPM_GRID

#if HAVE_DUNE_POLYGONGRID
#include <dune/polygongrid/declaration.hh>
#endif // #if HAVE_DUNE_POLYGONGRID

#include <dune/fem/misc/griddeclaration.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  struct HostGridAccess;



  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid >
    struct BoundaryIdProvider;



    // BoundaryIdProvider for AlbertaGrid
    // ----------------------------------

#if HAVE_ALBERTA
    template< int dim, int dimW >
    struct BoundaryIdProvider< AlbertaGrid< dim, dimW > >
    {
      typedef AlbertaGrid< dim, dimW > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_ALBERTA



    // BoundaryIdProvider for ALUGrid
    // ------------------------------

#if HAVE_DUNE_ALUGRID
    template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct BoundaryIdProvider< ALUGrid< dim, dimw, elType, refineType, Comm > >
    {
      typedef ALUGrid< dim, dimw, elType, refineType, Comm > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_DUNE_ALUGRID



    // BoundaryIdProvider for CacheItGrid
    // ----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< CacheItGrid< HostGrid > >
    {
      typedef CacheItGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for CartesianGrid
    // ------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< CartesianGrid< HostGrid > >
    {
      typedef CartesianGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for FilteredGrid
    // -----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< FilteredGrid< HostGrid > >
    {
      typedef FilteredGrid< HostGrid > GridType;

      // todo: FilteredGrid is a filtering grid and, hence, needs a specialized
      //       version of boundaryId.
      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        if( !HostGridAccess< GridType >::hostIntersection( intersection ).boundary() )
          DUNE_THROW( NotImplemented, "BoundaryIdProvider for artificial boundaries of FilteredGrid not implemented." );
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for GeometryGrid
    // -----------------------------------

    template< class HostGrid, class CoordFunction, class Allocator >
    struct BoundaryIdProvider< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Allocator > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };



    // BoundaryIdProvider for IdGrid
    // -----------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< IdGrid< HostGrid > >
    {
      typedef IdGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for OneDGrid
    // -------------------------------

    template<>
    struct BoundaryIdProvider< OneDGrid >
    {
      typedef OneDGrid GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.boundarySegmentIndex();
      }
    };



    // BoundaryIdProvider for ParallelGrid
    // -----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< ParallelGrid< HostGrid > >
    {
      typedef ParallelGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SphereGrid
    // ---------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid, class MapToSphere >
    struct BoundaryIdProvider< SphereGrid< HostGrid, MapToSphere > >
    {
      typedef SphereGrid< HostGrid, MapToSphere > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SPGrid
    // -----------------------------

#if HAVE_DUNE_SPGRID
    template< class ct, int dim, template< int > class Strategy, class Comm >
    struct BoundaryIdProvider< SPGrid< ct, dim, Strategy, Comm > >
    {
      typedef SPGrid< ct, dim, Strategy, Comm > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };
#endif // #if HAVE_DUNE_SPGRID


    // BoundaryIdProvider for PolygonGrid
    // ----------------------------------

#if HAVE_DUNE_POLYGONGRID
    template< class ct >
    struct BoundaryIdProvider< PolygonGrid< ct > >
    {
      typedef PolygonGrid< ct > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.impl().boundaryId()) : 0);
      }
    };
#endif // #if HAVE_OPM_GRID

    // BoundaryIdProvider for PolyhedralGrid
    // -------------------------------------

#if HAVE_OPM_GRID
    template< int dim, int dimworld, class ct >
    struct BoundaryIdProvider< PolyhedralGrid< dim, dimworld, ct > >
    {
      typedef PolyhedralGrid< dim, dimworld, ct > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.impl().boundaryId()) : 0);
      }
    };
#endif // #if HAVE_OPM_GRID


    // BoundaryIdProvider for UGGrid
    // -----------------------------

    template< int dim >
    struct BoundaryIdProvider< UGGrid< dim > >
    {
      typedef UGGrid< dim > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.boundarySegmentIndex();
      }
    };



    // BoundaryIdProvider for YaspGrid
    // -------------------------------

    template< int dim, class CoordCont >
    struct BoundaryIdProvider< YaspGrid< dim, CoordCont > >
    {
      typedef YaspGrid< dim, CoordCont > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };

    template< class GridPart, class Intersection>
    inline static int boundaryId ( const Intersection &intersection )
    {
      return Dune::Fem::BoundaryIdProvider< typename GridPart::GridType > ::
             boundaryId( intersection );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
