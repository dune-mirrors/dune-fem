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

    template< class GridPart, class Holder, class AliasType >
    void registerGrid ( pybind11::module module, pybind11::class_<GridPart,Holder,AliasType> &cls )
    {
      typedef typename GridPart::GridType HGrid;
      registerHierarchicalGrid< HGrid >( module );
      detail::registerGridPart< GridPart >( module, cls );
      cls.def( "coordinates", [] ( const GridPart &gridPart ) {
          return coordinates( static_cast< typename GridPart::GridViewType >( gridPart ) );
        } );
      cls.def( "tesselate", [] ( const GridPart &gridPart ) {
          return tesselate( static_cast< typename GridPart::GridViewType >( gridPart ) );
        } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
