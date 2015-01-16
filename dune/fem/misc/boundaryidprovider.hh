#ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
#define DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#error "Experimental grid extensions required. Reconfigure with --enable-experimental-grid-extensions."
#endif // #if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/declaration.hh>
#endif // #if HAVE_DUNE_SPGRID

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

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
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
#endif // #if HAVE_ALUGRID || HAVE_DUNE_ALUGRID



    // BoundaryIdProvider for ALUConformGrid
    // -------------------------------------

#if HAVE_ALUGRID
    template< int dim, int dimw >
    struct BoundaryIdProvider< ALUConformGrid< dim, dimw > >
    {
      typedef ALUConformGrid< dim, dimw > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_ALUGRID



    // BoundaryIdProvider for ALUCubeGrid
    // ----------------------------------

#if HAVE_ALUGRID
    template< int dim, int dimw >
    struct BoundaryIdProvider< ALUCubeGrid< dim, dimw > >
    {
      typedef ALUCubeGrid< dim, dimw > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_ALUGRID



    // BoundaryIdProvider for ALUSimplexGrid
    // -------------------------------------

#if HAVE_ALUGRID
    template< int dim, int dimw >
    struct BoundaryIdProvider< ALUSimplexGrid< dim, dimw > >
    {
      typedef ALUSimplexGrid< dim, dimw > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_ALUGRID



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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
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
        if( !iHostGridAccess< GridType >::getIntersection( intersection ).boundary() )
          DUNE_THROW( NotImplemented, "BoundaryIdProvider for artificial boundaries of FilteredGrid not implemented." );
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SGrid
    // ----------------------------

    template< int dim, int dimw, class ctype >
    struct BoundaryIdProvider< SGrid< dim, dimw, ctype > >
    {
      typedef SGrid< dim, dimw, ctype > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };



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
          ::boundaryId ( HostGridAccess< GridType >::getIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SPGrid
    // -----------------------------

#if HAVE_DUNE_SPGRID
    template< class ct, int dim, template< int > class Strategy, class Comm >
    struct BoundaryIdProvider< SPGrid< ct, dim, Strategy, Comm > >
    {
      typedef SPGrid< ct, dim, strategy, Comm > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };
#endif // #if HAVE_DUNE_SPGRID



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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
