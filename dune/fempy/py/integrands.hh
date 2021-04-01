#ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
#define DUNE_FEMPY_PY_INTEGRANDS_HH

#include <functional>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/python/common/typeregistry.hh>

#include <dune/fem/schemes/integrands.hh>

#include <dune/fempy/py/integrandsbase.hh>

#include <dune/fempy/pybind11/pybind11.hh>

// version 1.1

namespace Dune
{

  namespace FemPy
  {

    // registerIntegrands
    // ------------------

    template< class Integrands, class... options >
    inline void registerIntegrands ( pybind11::handle scope,
        pybind11::class_<Integrands,options...> cls)
    {
      typedef typename Integrands::GridPartType GridPart;
      typedef typename Integrands::DomainValueType DomainValue;
      typedef typename Integrands::RangeValueType RangeValue;
      typedef Fem::VirtualizedIntegrands< GridPart, DomainValue, RangeValue > VirtualizedIntegrands;

      detail::registerIntegrands( scope, cls );

      detail::clsVirtualizedIntegrands< VirtualizedIntegrands >( scope ).
        def( pybind11::init( [] ( Integrands &integrands ) {
          return new VirtualizedIntegrands( std::ref( integrands ) );
        }), pybind11::keep_alive< 1, 2 >() );
      pybind11::implicitly_convertible< Integrands, VirtualizedIntegrands >();
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
