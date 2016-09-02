#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <string>
#include <utility>

#include <dune/corepy/grid/gridview.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/hierarchical.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGrid
    // ------------

    template< class GridPart, class Holder, class AliasType >
    void registerGrid ( pybind11::module module, pybind11::class_<GridPart,Holder,AliasType> &cls )
    {
      // register a constructor and some properties that would not be carried over from
      // the GridView registration
      detail::registerGridPart< GridPart >( module, cls );
      registerHierarchicalGrid< typename GridPart::GridType >( module );
      auto clsGW = pybind11::class_<typename GridPart::GridViewType>(module, "GridView");
      clsGW.def(pybind11::init<GridPart&>());
      pybind11::implicitly_convertible<GridPart,typename GridPart::GridViewType>();
      CorePy::registerGridView<typename GridPart::GridViewType>( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
