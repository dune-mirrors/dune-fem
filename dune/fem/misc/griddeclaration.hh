#ifndef DUNE_FEM_GRIDDECLARATION_HH
#define DUNE_FEM_GRIDDECLARATION_HH

#include <dune/grid/geometrygrid/declaration.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/declaration.hh>
#endif // #if HAVE_DUNE_ALUGRID

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/declaration.hh>
#endif // #if HAVE_DUNE_SPGRID

#if HAVE_DUNE_METAGRID
#include <dune/grid/cacheitgrid/declaration.hh>
#include <dune/grid/cartesiangrid/declaration.hh>
#include <dune/grid/filteredgrid/declaration.hh>
#include <dune/grid/idgrid/declaration.hh>
#include <dune/grid/multispgrid/declaration.hh>
#include <dune/grid/parallelgrid/declaration.hh>
#include <dune/grid/spheregrid/declaration.hh>
#endif // #if HAVE_DUNE_METAGRID

namespace Dune
{

  // Forward Declarations for all Standard Grids
  // -------------------------------------------

  class OneDGrid;

  template< int dim >
  class UGGrid;

  template< int dim, class CoordCont >
  class YaspGrid;

  template< int dim, int dimworld >
  class AlbertaGrid;

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDDECLARATION_HH
