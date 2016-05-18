#ifndef DUNE_FEMPY_PYGRIDFUNCTION_HH
#define DUNE_FEMPY_PYGRIDFUNCTION_HH

#include <string>

#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pyvtk.hh>

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

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYGRIDFUNCTION_HH
