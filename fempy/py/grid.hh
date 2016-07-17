#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <string>
#include <utility>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/hierarchical.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGrid
    // ------------



    template< class GridPart >
    void registerGrid ( pybind11::module module )
    {
      typedef  typename GridPart::GridType HGrid;
      registerHierarchicalGrid< HGrid >( module );

      auto grid = registerGridPart< GridPart >( module, "LeafGrid" );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
