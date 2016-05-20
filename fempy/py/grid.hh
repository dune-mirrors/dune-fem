#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <string>
#include <utility>

#include <dune/fempy/py/grid/hierarchical.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGrid
    // ------------

    template< class GridPart >
    void registerGrid ( pybind11::module module )
    {
      static auto common = pybind11::module::import( "dune.common" );
      static auto femmpi = pybind11::module::import( "dune.femmpi" );

      registerHierarchicalGrid< HierarchicalGrid< typename GridPart::GridType > >( module );
      module.def( "makeSimplexGrid", &makeSimplexGrid< HierarchicalGrid< typename GridPart::GridType > > );

      registerGridPart< GridPart >( module, "LeafGrid" );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
