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

        // CoefficientsRegistry
        // --------------------

        template< class Integrands >
        struct CoefficientsRegistry
        {
          typedef std::array< pybind11::object, std::tuple_size< typename Integrands::CoefficientTupleType >::value > Objects;

          Objects &operator() ( pybind11::handle integrands )
          {
            auto result = objects_.insert( std::make_pair( integrands.ptr(), Objects() ) );
            auto pos = result.first;
            if( result.second )
            {
              pybind11::cpp_function remove_objects( [ this, pos ] ( pybind11::handle weakref ) {
                  objects_.erase( pos );
                  weakref.dec_ref();
                } );
              pybind11::weakref weakref( integrands, remove_objects );
              weakref.release();
            }
            return pos->second;
          }

          Objects &operator() ( Integrands &integrands )
          {
            return (*this)( pybind11::detail::get_object_handle( &integrands, pybind11::detail::get_type_info( typeid( Integrands ) ) ) );
          }

        private:
          std::map< PyObject *, Objects > objects_;
        };



        // coefficientsRegistry
        // --------------------

        template< class Integrands >
        DUNE_EXPORT inline CoefficientsRegistry< Integrands > &coefficientsRegistry ()
        {
          static CoefficientsRegistry< Integrands > registry;
          return registry;
        }




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



        // setCoefficient
        // --------------

        template< class Integrands, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::CoefficientTupleType >::value > 0) >
        setCoefficient ( pybind11::class_< Integrands, options... > cls, PriorityTag< 1 > )
        {
          const std::size_t numCoefficients = std::tuple_size< typename Integrands::CoefficientTupleType >::value;
          std::array< std::function< void( Integrands &, pybind11::handle ) >, numCoefficients > dispatch;
          Hybrid::forEach( std::make_index_sequence< numCoefficients >(), [ &dispatch ] ( auto i ) {
              dispatch[ i ] = [ i ] ( Integrands &integrands, pybind11::handle o ) {
                  integrands.template setCoefficient< i >( o.template cast< std::tuple_element_t< i, typename Integrands::CoefficientTupleType > >() );
                  coefficientsRegistry< Integrands >()( integrands )[ i ] = pybind11::object( o, true );
                };
            } );

          using pybind11::operator""_a;

          cls.def( "setCoefficient", [ dispatch ] ( Integrands &self, int k, pybind11::handle coefficient ) {
              if( k >= dispatch.size() )
                throw pybind11::index_error( "No such coefficient: " + std::to_string( k ) + " (must be in [0, " + std::to_string( dispatch.size() ) + "[)" );
              dispatch[ k ]( self, coefficient );
            }, "k"_a, "coefficient"_a );
        }

        template< class Integrands, class... options >
        inline static std::enable_if_t< (std::tuple_size< typename Integrands::CoefficientTupleType >::value == 0) >
        setCoefficient ( pybind11::class_< Integrands, options... > cls, PriorityTag< 1 > )
        {
          using pybind11::operator""_a;

          cls.def( "setCoefficient", [] ( Integrands &, int k, pybind11::handle ) {
              throw pybind11::index_error( "No such coefficient: " + std::to_string( k ) + " (there are no constants)" );
            }, "k"_a, "coefficient"_a );
        }

        template< class Integrands, class... options >
        inline static void setCoefficient ( pybind11::class_< Integrands, options... > cls, PriorityTag< 0 > )
        {}

      } // namespace RegisterIntegrands



      // registerIntegrands
      // ------------------

      template< class Integrands >
      inline void registerIntegrands ( pybind11::handle scope, pybind11::class_<Integrands> cls )
      {
        RegisterIntegrands::setConstant( cls, PriorityTag< 42 >() );
        RegisterIntegrands::setCoefficient( cls, PriorityTag< 42 >() );
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
