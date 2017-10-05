#ifndef DUNE_FEM_CAPABILITIES_HH
#define DUNE_FEM_CAPABILITIES_HH

#include <dune/common/version.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/misc/griddeclaration.hh>
#include <dune/fem/misc/metaprogramming.hh>

namespace Dune
{

  // Core Capabilities
  // -----------------

  namespace Capabilities
  {

    // hasHierarchicIndexSet
    // ---------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class Grid >
    struct hasHierarchicIndexSet< const Grid >
    {
      static const bool v = false;
    };

    template< class Grid >
    struct hasHierarchicIndexSet
    {
      static const bool v = false;
    };

#if HAVE_DUNE_ALUGRID
    template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct hasHierarchicIndexSet< ALUGrid< dim, dimw, elType, refineType, Comm > >
    {
      static const bool v = true;
    };
#endif // #if HAVE_DUNE_ALUGRID

    template<>
    struct hasHierarchicIndexSet< OneDGrid >
    {
      static const bool v = false;
    };

    template< int dim >
    struct hasHierarchicIndexSet< UGGrid< dim > >
    {
      static const bool v = false;
    };

    template< int dim, class CoordCont >
    struct hasHierarchicIndexSet< YaspGrid< dim, CoordCont > >
    {
      static const bool v = false;
    };


    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasHierarchicIndexSet< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      //static const bool v = hasHierarchicIndexSet< HostGrid > :: v;
      static const bool v = false ;
    };

  } // namespace Capabilities


  namespace Fem
  {

    // Fem Capabilities
    // ----------------

    namespace Capabilities
    {

      // hasAllCodimEntities
      // -------------------

      template< class Grid >
      class hasAllCodimEntities
      {
        template< unsigned int codim >
        struct Codim
        : public Dune::Capabilities::hasEntity< Grid, codim >
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
        static const bool v = Dune::Fem::Capabilities::supportsCallbackAdaptation< Grid > :: v;
      };

#if HAVE_DUNE_ALUGRID
      template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
      struct supportsCallbackAdaptation< ALUGrid< dim, dimworld, elType, refineType, Comm > >
      {
        static const bool v = true;
      };
#endif // #if HAVE_DUNE_ALUGRID

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

#if HAVE_DUNE_METAGRID
      template< class HostGrid >
      struct supportsCallbackAdaptation< CartesianGrid< HostGrid > >
      {
        static const bool v = supportsCallbackAdaptation< HostGrid > :: v;
      };
#endif // #if HAVE_DUNE_METAGRID



      // isLocallyAdaptive
      // -----------------

      template< class Grid >
      struct isLocallyAdaptive
      {
        static const bool v = false;
      };

      template< class Grid >
      struct isLocallyAdaptive< const Grid >
      {
        static const bool v = Dune::Fem::Capabilities::isLocallyAdaptive< Grid > :: v;
      };

#if HAVE_DUNE_ALUGRID
      template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
      struct isLocallyAdaptive< ALUGrid< dim, dimworld, elType, refineType, Comm > >
      {
        static const bool v = true;
      };
#endif //#f HAVE_DUNE_ALUGRID

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

#if HAVE_DUNE_METAGRID
      template< class HostGrid >
      struct isLocallyAdaptive< CartesianGrid< HostGrid > >
      {
        static const bool v = isLocallyAdaptive< HostGrid > :: v;
      };
#endif // #if HAVE_DUNE_METAGRID

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CAPABILITIES_HH
