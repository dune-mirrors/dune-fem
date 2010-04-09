#ifndef DUNE_FEM_CAPABILITIES_HH
#define DUNE_FEM_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/onedgrid.hh>

#include <dune/fem/misc/metaprogramming.hh>

namespace Dune
{

  // Forward Declarations for all Standard Grids
  // -------------------------------------------

  template< int dim, int dimw >
  class ALUConformGrid;

  template< int dim, int dimw >
  class ALUCubeGrid;

  template< int dim, int dimw >
  class ALUSimplexGrid;

  template< int dim, int dimw, class ctype >
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

  }

}

#endif // #ifndef DUNE_FEM_CAPABILITIES_HH
