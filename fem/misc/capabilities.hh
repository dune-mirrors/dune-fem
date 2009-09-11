#ifndef DUNE_FEM_CAPABILITIES_HH
#define DUNE_FEM_CAPABILITIES_HH

#include <dune/grid/onedgrid.hh>

namespace Dune
{

  // Forward Declarations for all Standard Grids
  // -------------------------------------------

  template< int dim, int dimw >
  class AlbertaGrid;

  template< int dim, int dimw >
  class ALUConformGrid;

  template< int dim, int dimw >
  class ALUCubeGrid;

  template< int dim, int dimw >
  class ALUSimplexGrid;

  template< int dim, int dimw >
  class SGrid;

  template< int dim >
  class UGGrid;

  template< int dim >
  class YaspGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class Grid >
    struct hasHierarchicIndexSet< const Grid >
    {
      static const bool v = hasHierarchicIndexSet< Grid >::v;
    };


    template< int dim, int dimw >
    class hasHierarchicIndexSet< AlbertaGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw >
    class hasHierarchicIndexSet< ALUConformGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw >
    class hasHierarchicIndexSet< ALUCubeGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template< int dim, int dimw >
    class hasHierarchicIndexSet< ALUSimplexGrid< dim, dimw > >
    {
      static const bool v = true;
    };

    template<>
    class hasHierarchicIndexSet< OneDGrid >
    {
      static const bool v = false;
    };

    template< int dim, int dimw >
    class hasHierarchicIndexSet< SGrid< dim, dimw > >
    {
      static const bool v = false;
    };

    template< int dim >
    class hasHierarchicIndexSet< UGGrid< dim > >
    {
      static const bool v = false;
    };

    template< int dim >
    class hasHierarchicIndexSet< YaspGrid< dim > >
    {
      static const bool v = false;
    };

  }

}

#endif // #ifndef DUNE_FEM_CAPABILITIES_HH
