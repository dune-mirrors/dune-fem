#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerSpace
    // -------------

    template< class Space >
    void registerSpace ( pybind11::module module )
    {
      typedef typename Space::GridPartType GridPart;

      pybind11::class_< Space > cls( module, "Space" );

      cls.def( "__init__", [] ( Space &instance, GridPart &grid ) {
          new( &instance ) Space( grid );
        }, pybind11::keep_alive< 1, 2 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
