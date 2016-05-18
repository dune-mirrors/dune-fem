#ifndef DUNE_FEMPY_PY_GRID_ENTITY_HH
#define DUNE_FEMPY_PY_GRID_ENTITY_HH

#include <string>
#include <tuple>
#include <utility>

#include <dune/fempy/py/grid/geometry.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGridEntity
    // ------------------

    template< class Entity >
    pybind11::class_< Entity > registerGridEntity ( pybind11::handle scope )
    {
      registerGridGeometry< typename Entity::Geometry >( scope );

      static const std::string clsName = "Entity" + std::to_string( Entity::codimension );
      pybind11::class_< Entity > cls( scope, clsName.c_str() );

      cls.def_property_readonly_static( "codimension", [] () { return  Entity::codimension; } );

      cls.def_property_readonly( "geometry", &Entity::geometry );
      cls.def_property_readonly( "level", &Entity::level );

      return cls;
    }



    // registerGridEntities
    // --------------------

    template< class Grid, int... codim >
    auto registerGridEntities ( pybind11::handle scope, std::integer_sequence< int, codim... > )
    {
      return std::make_tuple( registerGridEntity< typename Grid::template Codim< codim >::Entity >( scope )... );
    }

    template< class Grid >
    auto registerGridEntities ( pybind11::handle scope )
    {
      return registerGridEntities< Grid >( scope, std::make_integer_sequence< int, Grid::dimension+1 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_ENTITY_HH
