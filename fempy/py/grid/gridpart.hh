#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <string>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/misc/l2norm.hh>

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
      template <class GridPart, class Field, int dimR>
      double l2Norm (pybind11::object &o)
      {
        typedef Dune::FemPy::VirtualizedGridFunction<GridPart, Dune::FieldVector<Field,dimR>> VF;
        const VF &u = o.cast<VF>();
        Dune::Fem::L2Norm<GridPart> norm( u.gridPart(), 8);
        return norm.norm(u);
      }
      template< class GridPart, int... dimRange >
      auto defL2Norm ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
      {
        typedef std::function< double(pybind11::object&) > Dispatch;
        std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( l2Norm< GridPart, double, dimRange > )... }};

        return [ dispatch ] ( pybind11::object gp, pybind11::object func ) {
          const GridPart &gridPart = gp.cast< const GridPart & >();
          pybind11::object dimRobj = func.attr("dimRange");
          int dimR = func.attr("dimRange").cast<int>();
          if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ static_cast< std::size_t >( dimR ) ]( func );
        };
      }

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

        registerGridPartConstructorFromGrid<GridPart,Cls>(cls,std::is_constructible<GridPart,Grid&>());
        cls.attr( "dimGrid" ) = pybind11::int_( GridPart::dimension );
        cls.attr( "dimWorld" ) = pybind11::int_( GridPart::dimensionworld );

        cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
        cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

        cls.def( "l2Norm", defL2Norm< GridPart >( cls, "l2Norm", std::make_integer_sequence< int, 11 >() ) );
      }
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
