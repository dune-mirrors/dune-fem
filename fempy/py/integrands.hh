#ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
#define DUNE_FEMPY_PY_INTEGRANDS_HH

#include <functional>

#include <dune/fem/schemes/integrands.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerIntegrands
      // ------------------

      template< class Integrands >
      inline pybind11::class_< Integrands > registerIntegrands ( pybind11::handle scope, const char *clsName )
      {
        pybind11::class_< Integrands > cls( scope, clsName );
        return cls;
      }



      // clsVirtualizedIntegrands
      // ------------------------

      template< class GridPart, class Value >
      inline pybind11::class_< Fem::VirtualizedIntegrands< GridPart, Value > > clsVirtualizedIntegrands ( pybind11::handle scope )
      {
        typedef Fem::VirtualizedIntegrands< GridPart, Value > Integrands;
        static pybind11::class_< Integrands > cls = registerIntegrands< Integrands >( scope, "VirtualizedIntegrands" );
        return cls;
      }

    } // namespace detail



    // registerIntegrands
    // ------------------

    template< class Integrands >
    inline pybind11::class_< Integrands > registerIntegrands ( pybind11::handle scope, const char *clsName = "Integrands" )
    {
      typedef typename Integrands::GridPartType GridPart;
      typedef typename Integrands::ValueType Value;
      typedef Fem::VirtualizedIntegrands< GridPart, Value > VirtualizedIntegrands;

      pybind11::class_< Integrands > cls = detail::registerIntegrands< Integrands >( scope, clsName );

      detail::clsVirtualizedIntegrands< GridPart, Value >( scope ).def( "__init__", [] ( VirtualizedIntegrands &self, Integrands &integrands ) {
          new (&self) VirtualizedIntegrands( std::ref( integrands ) );
        }, pybind11::keep_alive< 1, 2 >() );
      pybind11::implicitly_convertible< Integrands, VirtualizedIntegrands >();

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
