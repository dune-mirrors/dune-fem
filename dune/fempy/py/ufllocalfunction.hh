#ifndef DUNE_FEMPY_PY_UFLLOCALFUNCTION_HH
#define DUNE_FEMPY_PY_UFLLOCALFUNCTION_HH

#include <functional>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/python/common/typeregistry.hh>

#include <dune/fempy/py/grid/function.hh>

#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      namespace RegisterUFLLocalFunction
      {

        // setConstant
        // -----------

        template< class UFLLocalFunction, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename UFLLocalFunction::ConstantTupleType >::value > 0) >
        setConstant ( pybind11::class_< UFLLocalFunction, options... > cls, PriorityTag< 1 > )
        {
          const std::size_t numConstants = std::tuple_size< typename UFLLocalFunction::ConstantTupleType >::value;
          std::array< std::function< void( UFLLocalFunction &, pybind11::handle ) >, numConstants > dispatch;
          Hybrid::forEach( std::make_index_sequence< numConstants >(), [ &dispatch ] ( auto i ) {
              dispatch[ i ] = [ i ] ( UFLLocalFunction &integrands, pybind11::handle o ) {
                  integrands.template constant< i.value >() = o.template cast< typename UFLLocalFunction::template ConstantType< i.value > >();
                };
            } );

          cls.def( "setConstant", [ dispatch ] ( UFLLocalFunction &self, int k, pybind11::handle o ) {
              if( static_cast< std::size_t >( k ) >= dispatch.size() )
                throw pybind11::index_error( "No such constant: " + std::to_string( k ) + " (must be in [0, " + std::to_string( dispatch.size() ) + "[)" );
              dispatch[ k ]( self, o );
            } );
        }

        template< class UFLLocalFunction, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename UFLLocalFunction::ConstantTupleType >::value == 0) >
        setConstant ( pybind11::class_< UFLLocalFunction, options... > cls, PriorityTag< 1 > )
        {
          cls.def( "setConstant", [] ( UFLLocalFunction &, int k, pybind11::handle ) {
              throw pybind11::index_error( "No such constant: " + std::to_string( k ) + " (there are no constants)" );
            } );
        }

        template< class UFLLocalFunction, class... options >
        inline static void setConstant ( pybind11::class_< UFLLocalFunction, options... > cls, PriorityTag< 0 > )
        {}

      } // namespace RegisterUFLLocalFunction



      // registerUFLLocalFunction
      // ------------------

      template< class UFLLocalFunction >
      inline void registerUFLLocalFunction ( pybind11::handle scope, pybind11::class_<UFLLocalFunction> cls )
      {
        RegisterUFLLocalFunction::setConstant( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerUFLLocalFunction
    // ------------------

    template< class UFLLocalFunction, class... options >
    inline void registerUFLLocalFunction ( pybind11::handle scope,
        pybind11::class_<UFLLocalFunction,options...> cls)
    {
      FemPy::registerGridFunction( scope, cls ); // could be confused with Dune::Python version
      detail::registerUFLLocalFunction( scope, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_UFLLOCALFUNCTIONS_HH
