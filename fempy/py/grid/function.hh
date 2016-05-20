#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <string>

#include <dune/fempy/py/grid/vtk.hh>
#include <dune/fempy/pybind11/pybind11.h>

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

        pybind11::class_< GridFunction > cls( scope, clsName );

        registerLocalFunction< LocalFunction >( cls );

        cls.def( "__repr__", [] ( GridFunction &gf ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + gf.name() + ")";
          } );

        cls.def_property_readonly( "dimRange", [] ( GridFunction &gf ) -> int { return GridFunction::RangeType::dimension; } );

        cls.def_property_readonly( "name", [] ( GridFunction &gf ) -> std::string { return gf.name(); } );

        cls.def( "localFunction", [] ( const GridFunction &gf, const Entity &entity ) -> LocalFunction {
            return gf.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >() );

        cls.def( "addToVTKWriter", &addToVTKWriter< GridFunction >, pybind11::keep_alive< 1, 2 >() );

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

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
