#ifndef DUNE_FEMPY_PY_INTEGRANDSBASE_HH
#define DUNE_FEMPY_PY_INTEGRANDSBASE_HH

#include <functional>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/python/common/typeregistry.hh>

#include <dune/fem/schemes/integrands.hh>

#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      namespace RegisterIntegrands
      {

        // setConstant
        // -----------

        template< class Integrands, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::ConstantTupleType >::value > 0) >
        setConstant ( pybind11::class_< Integrands, options... > cls, PriorityTag< 1 > )
        {
          const std::size_t numConstants = std::tuple_size< typename Integrands::ConstantTupleType >::value;
          std::array< std::function< void( Integrands &, pybind11::handle ) >, numConstants > dispatch;
          Hybrid::forEach( std::make_index_sequence< numConstants >(), [ &dispatch ] ( auto i ) {
              dispatch[ i ] = [ i ] ( Integrands &integrands, pybind11::handle o ) {
                  integrands.template constant< i.value >() = o.template cast< typename Integrands::template ConstantType< i.value > >();
                };
            } );

          cls.def( "setConstant", [ dispatch ] ( Integrands &self, int k, pybind11::handle o ) {
              if( static_cast< std::size_t >( k ) >= dispatch.size() )
                throw pybind11::index_error( "No such constant: " + std::to_string( k ) + " (must be in [0, " + std::to_string( dispatch.size() ) + "[)" );
              dispatch[ k ]( self, o );
            } );
        }

        template< class Integrands, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::ConstantTupleType >::value == 0) >
        setConstant ( pybind11::class_< Integrands, options... > cls, PriorityTag< 1 > )
        {
          cls.def( "setConstant", [] ( Integrands &, int k, pybind11::handle ) {
              throw pybind11::index_error( "No such constant: " + std::to_string( k ) + " (there are no constants)" );
            } );
        }

        template< class Integrands, class... options >
        inline static void setConstant ( pybind11::class_< Integrands, options... > cls, PriorityTag< 0 > )
        {}

      } // namespace RegisterIntegrands



      // registerIntegrands
      // ------------------

      template< class Integrands >
      inline void registerIntegrands ( pybind11::handle scope, pybind11::class_<Integrands> cls )
      {
        RegisterIntegrands::setConstant( cls, PriorityTag< 42 >() );
      }



      // clsVirtualizedIntegrands
      // ------------------------

      // template< class GridPart, class DomainValue, class RangeValue >
      template< class VIntegrands >
      inline pybind11::class_< VIntegrands >
      clsVirtualizedIntegrands ( pybind11::handle scope )
      {
        auto cls = Python::insertClass<VIntegrands>(scope,"VirtualizedIntegrands",
            Python::GenerateTypeName("TODO") );
        if (cls.second)
          registerIntegrands( scope, cls.first );
        return cls.first;
      }

    } // namespace detail

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_INTEGRANDSBASE_HH
