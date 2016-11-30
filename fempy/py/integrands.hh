#ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
#define DUNE_FEMPY_PY_INTEGRANDS_HH

#include <functional>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typeutilities.hh>

#include <dune/fem/schemes/integrands.hh>

#include <dune/corepy/pybind11/pybind11.h>

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

        template< class Integrands, class Holder, class Alias >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::ConstantTupleType >::value > 0) >
        setConstant ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 1 > )
        {
          const std::size_t numConstants = std::tuple_size< typename Integrands::ConstantTupleType >::value;
          std::array< std::function< void( Integrands &, pybind11::handle ) >, numConstants > dispatch;
          Hybrid::forEach( std::make_index_sequence< numConstants >(), [ &dispatch ] ( auto i ) {
              dispatch[ i ] = [ i ] ( Integrands &integrands, pybind11::handle o ) {
                  integrands.template constant< i >() = o.template cast< typename Integrands::template ConstantType< i > >();
                };
            } );

          cls.def( "setConstant", [ dispatch ] ( Integrands &integrands, int k, pybind11::handle o ) {
              if( k >= dispatch.size() )
                throw std::range_error( "No such constant: " + std::to_string( k ) + " (must be in [0, " + std::to_string( dispatch.size() ) + "[)" );
              dispatch[ k ]( integrands, o );
            } );
        }

        template< class Integrands, class Holder, class Alias >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::ConstantTupleType >::value == 0) >
        setConstant ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 1 > )
        {
          cls.def( "setConstant", [] ( Integrands &integrands, int k, pybind11::handle ) {
              throw std::range_error( "No such constant: " + std::to_string( k ) + " (there are no constants)" );
            } );
        }

        template< class Integrands, class Holder, class Alias >
        inline static void setConstant ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 0 > )
        {}



        // setCoefficient
        // --------------

        template< class Integrands, class Holder, class Alias >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::CoefficientTupleType >::value > 0) >
        setCoefficient ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 1 > )
        {
          const std::size_t numCoefficients = std::tuple_size< typename Integrands::CoefficientTupleType >::value;
          std::array< std::function< void( Integrands &, pybind11::handle ) >, numCoefficients > dispatch;
          Hybrid::forEach( std::make_index_sequence< numCoefficients >(), [ &dispatch ] ( auto i ) {
              dispatch[ i ] = [ i ] ( Integrands &integrands, pybind11::handle o ) {
                  integrands.template setCoefficient< i >( o.template cast< std::tuple_element_t< i, typename Integrands::CoefficientTupleType > >() );
                };
            } );

          cls.def( "setCoefficient", [ dispatch ] ( Integrands &integrands, int k, pybind11::handle o ) {
              if( k >= dispatch.size() )
                throw std::range_error( "No such coefficient: " + std::to_string( k ) + " (must be in [0, " + std::to_string( dispatch.size() ) + "[)" );
              dispatch[ k ]( integrands, o );
            } );
        }

        template< class Integrands, class Holder, class Alias >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::CoefficientTupleType >::value == 0) >
        setCoefficient ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 1 > )
        {
          cls.def( "setCoefficient", [] ( Integrands &integrands, int k, pybind11::handle ) {
              throw std::range_error( "No such coefficient: " + std::to_string( k ) + " (there are no constants)" );
            } );
        }

        template< class Integrands, class Holder, class Alias >
        inline static void setCoefficient ( pybind11::class_< Integrands, Holder, Alias > cls, PriorityTag< 0 > )
        {}

      } // namespace RegisterIntegrands



      // registerIntegrands
      // ------------------

      template< class Integrands >
      inline pybind11::class_< Integrands > registerIntegrands ( pybind11::handle scope, const char *clsName )
      {
        pybind11::class_< Integrands > cls( scope, clsName );
        RegisterIntegrands::setConstant( cls, PriorityTag< 42 >() );
        RegisterIntegrands::setCoefficient( cls, PriorityTag< 42 >() );
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
