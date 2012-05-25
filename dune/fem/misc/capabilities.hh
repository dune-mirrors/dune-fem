#ifndef DUNE_FEM_CAPABILITIES_HH
#define DUNE_FEM_CAPABILITIES_HH

#include <dune/common/version.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/misc/metaprogramming.hh>

// include forward declarations of standard grids 
#include <dune/fem/misc/griddeclaration.hh>

namespace Dune
{
  // Capabilities
  // ------------

  namespace Capabilities
  {

    // hasHierarchicIndexSet
    // ---------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class Grid >
    struct hasHierarchicIndexSet< const Grid >
    {
      static const bool v = hasHierarchicIndexSet< Grid >::v;
    };


    template< int dim, int dimw >
    struct hasHierarchicIndexSet< ALUConformGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw >
    struct hasHierarchicIndexSet< ALUCubeGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw >
    struct hasHierarchicIndexSet< ALUSimplexGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct hasHierarchicIndexSet< ALUGrid< dim, dimw, elType, refineType, Comm > >
    {
      static const bool v = true;
    };

    template<>
    struct hasHierarchicIndexSet< OneDGrid >
    {
      static const bool v = false;
    };

    template< int dim, int dimw, class ctype >
    struct hasHierarchicIndexSet< SGrid< dim, dimw, ctype > >
    {
      static const bool v = false;
    };

    template< int dim >
    struct hasHierarchicIndexSet< UGGrid< dim > >
    {
      static const bool v = false;
    };

    template< int dim >
    struct hasHierarchicIndexSet< YaspGrid< dim > >
    {
      static const bool v = false;
    };


    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasHierarchicIndexSet< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      //static const bool v = hasHierarchicIndexSet< HostGrid > :: v;
      static const bool v = false ;
    };

    // hasAllCodimEntities
    // -------------------

    template< class Grid >
    class hasAllCodimEntities
    {
      template< unsigned int codim >
      struct Codim
      : public hasEntity< Grid, codim >
      {};
    
    public:
      static const bool v = Loop< MetaAnd, Codim, Grid :: dimension > :: v;
      static const bool value = v;
    };



    // supportsCallbackAdaptation
    // --------------------------

    template< class Grid >
    struct supportsCallbackAdaptation
    {
      static const bool v = false;
    };

    template< class Grid >
    struct supportsCallbackAdaptation< const Grid >
    {
      static const bool v = Dune :: Capabilities :: supportsCallbackAdaptation< Grid > :: v;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUSimplexGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUCubeGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUConformGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct supportsCallbackAdaptation< ALUGrid< dim, dimworld, elType, refineType, Comm > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct supportsCallbackAdaptation< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = supportsCallbackAdaptation< HostGrid > :: v;
    };

    template< class HostGrid >
    struct supportsCallbackAdaptation< CartesianGrid< HostGrid > >
    {
      static const bool v = supportsCallbackAdaptation< HostGrid > :: v;
    };

    // isLocallyAdaptive
    // --------------------------

    template< class Grid >
    struct isLocallyAdaptive
    {
      static const bool v = false;
    };

    template< class Grid >
    struct isLocallyAdaptive< const Grid >
    {
      static const bool v = Dune :: Capabilities :: isLocallyAdaptive< Grid > :: v;
    };

    template< int dim, int dimworld >
    struct isLocallyAdaptive< ALUSimplexGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct isLocallyAdaptive< ALUCubeGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct isLocallyAdaptive< ALUConformGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct isLocallyAdaptive< ALUGrid< dim, dimworld, elType, refineType, Comm > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct isLocallyAdaptive< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim >
    struct isLocallyAdaptive< UGGrid< dim > >
    {
      static const bool v = true;
    };

    template<>
    struct isLocallyAdaptive< OneDGrid >
    {
      static const bool v = true;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct isLocallyAdaptive< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = isLocallyAdaptive< HostGrid > :: v;
    };

    template< class HostGrid >
    struct isLocallyAdaptive< CartesianGrid< HostGrid > >
    {
      static const bool v = isLocallyAdaptive< HostGrid > :: v;
    };


    // IsUnstructured (deprecated)
    // --------------
#ifdef DUNE_FEM_COMPATIBILITY
    template< class Grid >
    struct IsUnstructured
    {
      static const bool v = ! isCartesian< Grid > :: v ;
    };

    template< class Grid >
    struct IsUnstructured< const Grid >
    {
      static const bool v = ! isCartesian< Grid > :: v;
    };
#endif

  }

}

#endif // #ifndef DUNE_FEM_CAPABILITIES_HH
