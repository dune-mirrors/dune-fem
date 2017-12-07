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
                  integrands.template constant< i >() = o.template cast< typename Integrands::template ConstantType< i > >();
                };
            } );

          cls.def( "setConstant", [ dispatch ] ( Integrands &self, int k, pybind11::handle o ) {
              if( k >= dispatch.size() )
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

      template< class GridPart, class DomainValue, class RangeValue >
      inline pybind11::class_< Fem::VirtualizedIntegrands< GridPart, DomainValue, RangeValue > >
      clsVirtualizedIntegrands ( pybind11::handle scope )
      {
        typedef Fem::VirtualizedIntegrands< GridPart, DomainValue, RangeValue > Integrands;
        auto cls = Python::insertClass<Integrands>(scope,"VirtualizedIntegrads",
            Python::GenerateTypeName("TODO") );
        if (cls.second)
          registerIntegrands( scope, cls.first );
        return cls.first;
      }

    } // namespace detail



    // registerIntegrands
    // ------------------

    template< class Integrands >
    inline pybind11::class_< Integrands > registerIntegrands ( pybind11::handle scope, const char *clsName = "Integrands" )
    {
      typedef typename Integrands::GridPartType GridPart;
      typedef typename Integrands::DomainValueType DomainValue;
      typedef typename Integrands::RangeValueType RangeValue;
      typedef Fem::VirtualizedIntegrands< GridPart, DomainValue, RangeValue > VirtualizedIntegrands;

      auto cls = Python::insertClass<Integrands>(scope, clsName,
          Python::GenerateTypeName("TODO") );
      if (cls.second)
        detail::registerIntegrands( scope, cls.first );

      detail::clsVirtualizedIntegrands< GridPart, DomainValue, RangeValue >( scope ).
        def( pybind11::init( [] ( Integrands &integrands ) {
          return new VirtualizedIntegrands( std::ref( integrands ) );
        }), pybind11::keep_alive< 1, 2 >() );
      pybind11::implicitly_convertible< Integrands, VirtualizedIntegrands >();

      return cls.first;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_INTEGRANDS_HH
