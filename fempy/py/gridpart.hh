#ifndef DUNE_FEMPY_PY_GRIDPART_HH
#define DUNE_FEMPY_PY_GRIDPART_HH

#include <string>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/range.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGridPart
    // ----------------

    template< class Holder, class AliasType, class ...Args >
    void registerGridPart ( pybind11::handle scope, pybind11::class_<Dune::Fem::GeometryGridPart<Args...>,Holder,AliasType> &cls )
    {
      typedef Dune::Fem::GeometryGridPart<Args...> GridPart;
      typedef typename GridPart::GridFunctionType GridFunctionType;

      detail::registerGridPart<GridPart>(scope,cls);
      cls.def( "__init__", [] ( GridPart &instance,  GridFunctionType &gf ) {
         new (&instance) GridPart( gf );
      }, pybind11::keep_alive< 1, 2 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
