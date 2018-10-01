#ifndef DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH
#define DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH

#include <dune/fem/gridpart/common/gridpartadapter.hh>

namespace Dune
{

  namespace FemPy
  {

    template < class GridView >
    using GridPartAdapter = Fem::GridPartAdapter< GridView >;

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_GRIDPARTADAPTER_HH
