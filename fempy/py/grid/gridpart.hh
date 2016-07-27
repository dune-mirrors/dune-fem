#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <string>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/numpy.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerGridPart
    // ----------------

    namespace detail
    {
      template<class GridPart, class Cls>
      void registerGridPartConstructorFromGrid(Cls &cls, std::false_type)
      {}
      template<class GridPart, class Cls>
      void registerGridPartConstructorFromGrid(Cls &cls, std::true_type)
      {
        typedef typename GridPart::GridType Grid;
        cls.def("__init__",
            [] (GridPart &instance, Grid &grid) {
              new (&instance) GridPart(grid);
            },
            pybind11::keep_alive<1, 2>());
      }
      template< class GridPart, class Cls >
      void registerGridPart ( pybind11::handle scope, Cls &cls )
      {
        typedef typename GridPart::GridType Grid;

        const int dim = GridPart::dimension;

        registerGridPartConstructorFromGrid<GridPart,Cls>(cls,std::is_constructible<GridPart,Grid&>());
        cls.attr( "dimGrid" ) = pybind11::int_( GridPart::dimension );
        cls.attr( "dimWorld" ) = pybind11::int_( GridPart::dimensionworld );

        cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
        cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );
      }
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
