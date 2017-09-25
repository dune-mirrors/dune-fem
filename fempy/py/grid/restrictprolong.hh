#ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
#define DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH

#include <dune/common/visibility.hh>

#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // registerRestrictProlong
    // -----------------------

    template< class RestrictProlong >
    inline static pybind11::class_< RestrictProlong > registerRestrictProlong ( pybind11::handle scope, const char *name = "RestrictProlong" )
    {
      pybind11::class_< RestrictProlong > cls( scope, name );
      return cls;
    }


    namespace detail
    {

      // clsVirtualizedRestrictProlong
      // -----------------------------

      template< class Grid >
      DUNE_EXPORT inline pybind11::class_< VirtualizedRestrictProlong< Grid > > clsVirtualizedRestrictProlong ( pybind11::handle scope )
      {
        typedef VirtualizedRestrictProlong< Grid > RestrictProlong;
        static pybind11::class_< RestrictProlong > cls = registerRestrictProlong< RestrictProlong >( scope );
        return cls;
      }

    } // namespace detail

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
