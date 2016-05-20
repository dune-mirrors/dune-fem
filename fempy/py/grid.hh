#ifndef DUNE_FEMPY_PY_GRID_HH
#define DUNE_FEMPY_PY_GRID_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/hierarchical.hh>
#include <dune/fempy/py/grid/range.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/vtk.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
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



    // registerPyGlobalGridFunction
    // ----------------------------

    template< class GridPart, int dimRange >
    auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
    {
      typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
      static const std::string clsName = name + std::to_string( dimRange );
      return registerGridFunction< GridFunction >( scope, clsName.c_str() );
    };



    // registerPyLocalGridFunction
    // ---------------------------

    template< class GridPart, int dimRange >
    auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
    {
      typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
      static const std::string clsName = name + std::to_string( dimRange );
      return registerGridFunction< GridFunction >( scope, clsName.c_str() );
    };



    // defGlobalGridFunction
    // ---------------------

    template< class GridPart, int... dimRange >
    auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( registerPyGlobalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      auto makeDispatch = [] ( auto dimR ) {
          return [ dimR ] ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent ) {
              return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), std::move( evaluate ), dimR ), pybind11::return_value_policy::move, parent );
            };
        };
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( makeDispatch( std::integral_constant< int, dimRange >() ) )... }};

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



    // defLocalGridFunction
    // --------------------

    template< class GridPart, int... dimRange >
    auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( registerPyLocalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, pybind11::function, pybind11::object ) > Dispatch;
      auto makeDispatch = [] ( auto dimR ) {
          return [ dimR ] ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent ) {
              return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), std::move( evaluate ), dimR ), pybind11::return_value_policy::move, parent );
            };
        };
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( makeDispatch( std::integral_constant< int, dimRange >() ) )... }};

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



    // registerGrid
    // ------------

    template< class GridPart >
    pybind11::class_< GridPart > registerGrid ( pybind11::module module )
    {
      static auto common = pybind11::module::import( "dune.common" );
      static auto femmpi = pybind11::module::import( "dune.femmpi" );

      registerHierarchicalGrid< HierarchicalGrid< typename GridPart::GridType > >( module );
      module.def( "makeSimplexGrid", &makeSimplexGrid< HierarchicalGrid< typename GridPart::GridType > > );

      const int dim = GridPart::dimension;

      pybind11::class_< GridPart > cls( module, "LeafGrid" );
      cls.def( "__init__", [] ( GridPart &instance, HierarchicalGrid< typename GridPart::GridType > &hGrid ) {
          new (&instance) GridPart( *hGrid.grid() );
        }, pybind11::keep_alive< 1, 2 >() );

      registerPyGridPartRange< GridPart, 0 >( cls, "Elements" );
      cls.def_property_readonly( "elements", [] ( pybind11::object gridPart ) {
          return PyGridPartRange< GridPart, 0 >( gridPart.template cast< const GridPart & >(), gridPart );
        } );

      registerPyGridPartRange< GridPart, dim >( cls, "Vertices" );
      cls.def_property_readonly( "vertices", [] ( pybind11::object gridPart ) {
          return PyGridPartRange< GridPart, dim >( gridPart.template cast< const GridPart & >(), gridPart );
        } );

      cls.def( "__repr__", [] ( const GridPart &gridPart ) -> std::string {
          return "LeafGrid with " + std::to_string( gridPart.indexSet().size( 0 ) ) + " elements";
        } );

      cls.def_property_readonly( "hierarchicalGrid", [] ( GridPart &gridPart ) { return hierarchicalGrid( gridPart.grid() ); } );

      cls.def( "size", [] ( const GridPart &gridPart, int codim ) { return gridPart.indexSet().size( codim ); } );

      registerVTKWriter< typename GridPart::GridViewType >( module );
      cls.def( "vtkWriter", [] ( const GridPart &gridPart ) {
          return new VTKWriter< typename GridPart::GridViewType >( static_cast< typename GridPart::GridViewType >( gridPart ) );
        }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
      cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HH
