#ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
#define DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH

#include <dune/common/visibility.hh>

#include <dune/python/common/typeregistry.hh>

#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // clsVirtualizedRestrictProlong
      // -----------------------------

      template< class Grid >
      inline pybind11::class_< VirtualizedRestrictProlong< Grid > > clsVirtualizedRestrictProlong ( pybind11::handle scope )
      {
        typedef VirtualizedRestrictProlong< Grid > RestrictProlong;
        auto cls = Python::insertClass< RestrictProlong >(scope, "RestrictProlong",
            Python::GenerateTypeName("VirtualizedRestrictProlong",MetaType<Grid>()),
            Python::IncludeFiles{"dune/fempy/grid/virtualizedrestrictprolong.hh"});
        return cls.first;
      }

    } // namespace detail

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_RESTRICTPROLONG_HH
