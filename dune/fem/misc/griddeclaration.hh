#ifndef DUNE_FEM_GRIDDECLARATION
#define DUNE_FEM_GRIDDECLARATION

// include ALUGrid forward declaration 
#include <dune/grid/alugrid/common/declaration.hh>

namespace Dune {

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
#endif
