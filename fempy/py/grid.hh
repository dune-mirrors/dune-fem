#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <string>
#include <utility>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/hierarchical.hh>
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
      typedef  typename GridPart::GridType HGrid;
      registerHierarchicalGrid< HGrid >( module );

      module.def( "reader", [](const std::tuple<Reader,std::string> &constructor)
          { return reader<HGrid>(constructor); } );
      module.def( "reader", [](const std::string &constructor)
          { return reader<HGrid>(std::make_tuple(Reader::dgf,constructor)); } );
      module.def( "reader", [](const pybind11::dict &constructor)
          { return reader<HGrid>(constructor); } );


      auto grid = registerGridPart< GridPart >( module, "LeafGrid" );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
