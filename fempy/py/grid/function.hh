#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/numpy.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction >
    pybind11::class_< LocalFunction > registerLocalFunction ( pybind11::handle scope, const char *clsName = "LocalFunction" )
    {
      typedef typename LocalFunction::LocalCoordinateType LocalCoordinate;

      pybind11::class_< LocalFunction > cls( scope, clsName );

      cls.def_property_readonly( "dimRange", [] ( LocalFunction & ) -> int { return LocalFunction::RangeType::dimension; } );
      cls.def( "evaluate", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          lf.evaluate( x, value );
          return value;
        } );
      cls.def( "jacobian", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::JacobianRangeType jacobian;
          lf.jacobian( x, jacobian );
          return jacobian;
        } );

      return cls;
    }



    namespace detail
    {

      // registerGridFunction
      // --------------------

      template< class GridFunction >
      pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
      {
        typedef typename GridFunction::LocalFunctionType LocalFunction;
        typedef typename LocalFunction::EntityType Entity;
        typedef typename GridFunction::GridPartType GridPartType;

        pybind11::class_< GridFunction > cls( scope, clsName );

        registerLocalFunction< LocalFunction >( cls );

        cls.def( "__repr__", [] ( GridFunction &gf ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + gf.name() + ")";
          } );

        cls.def_property_readonly( "dimRange", [] ( GridFunction &gf ) -> int { return GridFunction::RangeType::dimension; } );

        cls.def_property_readonly( "name", [] ( GridFunction &gf ) -> std::string { return gf.name(); } );
        cls.def_property_readonly( "grid", [](GridFunction &gf) -> const GridPartType& {return gf.gridPart();} );

        cls.def( "localFunction", [] ( const GridFunction &gf, const Entity &entity ) -> LocalFunction {
            return gf.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

        cls.def( "addToVTKWriter", &Dune::CorePy::addToVTKWriter< GridFunction, VTKWriter< typename GridFunction::GridPartType::GridViewType > >,
            pybind11::keep_alive< 2, 1 >() );

        cls.def( "cellData", [] ( const GridFunction &gf ) { return cellData( gf ); } );
        cls.def( "pointData", [] ( const GridFunction &gf ) { return pointData( gf ); } );

        return cls;
      }



      // clsVirtualizedGridFunction
      // --------------------------

      template< class GridPart, class Value >
      inline pybind11::class_< VirtualizedGridFunction< GridPart, Value > > clsVirtualizedGridFunction ( pybind11::handle scope )
      {
        typedef VirtualizedGridFunction< GridPart, Value > GridFunction;
        static const std::string clsName = "VirtualizedGridFunction" + std::to_string( Value::dimension );
        static pybind11::class_< GridFunction > cls = registerGridFunction< GridFunction >( scope, clsName.c_str() );
        return cls;
      }

    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction >
    pybind11::class_< GridFunction > registerGridFunction ( pybind11::handle scope, const char *clsName = "GridFunction" )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Value;

      auto cls = detail::registerGridFunction< GridFunction >( scope, clsName );

      detail::clsVirtualizedGridFunction< GridPart, Value >( scope ).def( pybind11::init< GridFunction >() );
      pybind11::implicitly_convertible< GridFunction, VirtualizedGridFunction< GridPart, Value > >();

      return cls;
    }



    // registerVirtualizedGridFunction
    // -------------------------------

    template< class GridPart, int... dimRange >
    void registerVirtualizedGridFunction ( pybind11::handle scope, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::clsVirtualizedGridFunction< GridPart, FieldVector< double, dimRange > >( scope )... );
    };



    namespace detail
    {

      // makePyGlobalGridFunction
      // ------------------------

      template< class GridPart, int dimRange >
      auto makePyGlobalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( x ) );
            return v.template cast< FieldVector< double, dimRange > >();
          } );
      }



      // registerPyGlobalGridFunction
      // ----------------------------

      template< class GridPart, int dimRange >
      auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      pybind11::object pyGlobalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defGlobalGridFunction
    // ---------------------

    template< class GridPart, int... dimRange >
    auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyGlobalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyGlobalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gridPart, std::string name, pybind11::function evaluate ) {
          typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate x( 0 );
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate( x ) );
          const std::size_t dimR = len( v );
          if( dimR >= dispatch.size() )
            DUNE_THROW( NotImplemented, "globalGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ dimR ]( gridPart.cast< const GridPart & >(), std::move( name ), std::move( evaluate ), gridPart );
        };
    }



    namespace detail
    {

      // makePyLocalGridFunction
      // -----------------------

      template< class GridPart, int dimRange >
      auto makePyLocalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::EntityType Entity;
        typedef typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Entity &entity, const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( entity, x ) );
            return v.template cast< FieldVector< double, dimRange > >();
          } );
      }


      // registerPyLocalGridFunction
      // ---------------------------

      template< class GridPart, int dimRange >
      auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        static const std::string clsName = name + std::to_string( dimRange );
        return FemPy::registerGridFunction< GridFunction >( scope, clsName.c_str() );
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      pybind11::object pyLocalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defLocalGridFunction
    // --------------------

    template< class GridPart, int... dimRange >
    auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyLocalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyLocalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gp, std::string name, pybind11::function evaluate ) {
          const GridPart &gridPart = gp.cast< const GridPart & >();
          int dimR = -1;
          if( gridPart.template begin< 0 >() != gridPart.template end< 0 >() )
          {
            typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate x( 0 );
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( *gridPart.template begin< 0 >(), x ) );
            dimR = len( v );
          }
          dimR = gridPart.comm().max( dimR );
          if( dimR < 0 )
            DUNE_THROW( InvalidStateException, "Cannot create local grid function on empty grid" );
          if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ static_cast< std::size_t >( dimR ) ]( gridPart, std::move( name ), std::move( evaluate ), std::move( gp ) );
        };
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
