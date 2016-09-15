#ifndef DUNE_FEMPY_PY_GRIDPART_HH
#define DUNE_FEMPY_PY_GRIDPART_HH

#include <dune/corepy/common/typeregistry.hh>
#include <dune/corepy/grid/gridview.hh>

#include <dune/fempy/py/grid/gridpart.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGridPart
    // ----------------

    template< class GridPart, class Holder, class Alias >
    void registerGridPart ( pybind11::handle module, pybind11::class_< GridPart, Holder, Alias > &cls )
    {
      detail::registerGridPart< GridPart >( module, cls );

      typedef typename GridPart::GridViewType GridView;

      auto &registry = CorePy::typeRegistry();
      auto entry = registry.insertInner< GridPart, GridView >( "GridViewType", {} );
      registry.exportToPython( cls, entry.first->second );

      auto clsGW = pybind11::class_< GridView >( module, "GridView" );
      clsGW.def( pybind11::init< GridPart & >() );
      pybind11::implicitly_convertible< GridPart, GridView >();
      CorePy::registerGridView< GridView >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
