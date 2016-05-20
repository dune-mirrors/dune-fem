#ifndef DUNE_FEMPY_PY_FUNCTION_DISCRETE_HH
#define DUNE_FEMPY_PY_FUNCTION_DISCRETE_HH

#include <string>
#include <utility>

#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerDiscreteFunction
    // ------------------------

    template< class DiscreteFunction >
    pybind11::class_< DiscreteFunction > registerDiscreteFunction ( pybind11::handle scope, const char *clsName = "DiscreteFunction" )
    {
      auto cls = detail::registerGridFunction< DiscreteFunction >( scope, clsName );

      cls.def_property_readonly( "space", &DiscreteFunction::space );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_FUNCTION_DISCRETE_HH
