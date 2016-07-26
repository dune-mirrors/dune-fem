#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <string>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/range.hh>
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

        cls.attr( "dimGrid" ) = pybind11::int_( GridPart::dimension );
        cls.attr( "dimWorld" ) = pybind11::int_( GridPart::dimensionworld );

        registerGridPartConstructorFromGrid<GridPart,Cls>(cls,std::is_constructible<GridPart,Grid&>());

        registerPyGridPartRange< GridPart, 0 >( cls, "Elements" );
        cls.def_property_readonly( "elements", [] ( pybind11::object gridPart ) {
            return PyGridPartRange< GridPart, 0 >( gridPart.template cast< const GridPart & >(), gridPart );
          } );

        registerPyGridPartRange< GridPart, dim >( cls, "Vertices" );
        cls.def_property_readonly( "vertices", [] ( pybind11::object gridPart ) {
            return PyGridPartRange< GridPart, dim >( gridPart.template cast< const GridPart & >(), gridPart );
          } );

        cls.def( "__repr__", [] ( const GridPart &gridPart ) -> std::string {
            return "GridPart with " + std::to_string( gridPart.indexSet().size( 0 ) ) + " elements";
          } );

        cls.def_property_readonly( "hierarchicalGrid", [] ( GridPart &gridPart ) -> Grid & { return gridPart.grid(); } );

        cls.def( "size", [] ( const GridPart &gridPart, int codim ) { return gridPart.indexSet().size( codim ); } );

        /*
        cls.def( "coordinates", [] ( const GridPart &gridPart ) {
            return coordinates( static_cast< typename GridPart::GridViewType >( gridPart ) );
          } );
        cls.def( "tesselate", [] ( const GridPart &gridPart ) {
            return tesselate( static_cast< typename GridPart::GridViewType >( gridPart ) );
          } );
        */

        Dune::CorePy::registerVTKWriter< typename GridPart::GridViewType >( cls );

        cls.def( "vtkWriter", [] ( const GridPart &gridPart ) {
            return new VTKWriter< typename GridPart::GridViewType >( static_cast< typename GridPart::GridViewType >( gridPart ) );
          }, pybind11::keep_alive< 0, 1 >() );
        cls.def( "vtkWriter", [] ( const GridPart &gridPart, int level ) {
            return new SubsamplingVTKWriter< typename GridPart::GridViewType >( static_cast< typename GridPart::GridViewType >( gridPart ), level );
          }, pybind11::keep_alive< 0, 1 >() );

        cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
        cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

      }
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
