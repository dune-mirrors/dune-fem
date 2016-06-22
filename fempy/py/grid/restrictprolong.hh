#ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
#define DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerRestrictProlong
    // -----------------------

    template< class RestrictProlong >
    pybind11::class_< RestrictProlong > registerRestrictProlong ( pybind11::handle scope )
    {
      static pybind11::class_< RestrictProlong > cls( scope, "RestrictProlong" );
      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
