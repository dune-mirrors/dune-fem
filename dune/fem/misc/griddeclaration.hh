#ifndef DUNE_FEM_GRIDDECLARATION_HH
#define DUNE_FEM_GRIDDECLARATION_HH

// include ALUGrid forward declaration 
#if HAVE_DUNE_ALUGRID 
#include <dune/alugrid/common/declaration.hh>
#else
#include <dune/grid/alugrid/common/declaration.hh>
#endif

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

  class OneDGrid;

  template< int dim, int dimw, class ctype >
  class SGrid;

  template< int dim >
  class UGGrid;

  template< int dim >
  class YaspGrid;

  template< int dim, int dimworld >
  class AlbertaGrid;

  template< class HostGrid, class CoordFunction, class Allocator >
  class GeometryGrid;

  template< class HostGrid >
  class CartesianGrid;

}
#endif // #ifndef DUNE_FEM_GRIDDECLARATION_HH
