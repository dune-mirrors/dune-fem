#ifndef DUNE_FEMPY_PY_GRIDVIEW
#define DUNE_FEMPY_PY_GRIDVIEW
#include <dune/python/grid/gridview.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>

#include <dune/fem/marking/default.hh>
#include <dune/fem/marking/doerfler.hh>

namespace Dune
{
  namespace FemPy
  {
    template< class GridView, class... options >
    inline static void registerGridView ( pybind11::handle scope, pybind11::class_< GridView, options... > cls )
    {
      Dune::Python::registerGridView(scope,cls);
      typedef GridPart< GridView > GridPart;
      cls.def_property_readonly("canAdapt",[](GridView &self){
       return Dune::Fem::Capabilities::isAdaptiveIndexSet<typename GridPart::IndexSetType>::v;
          });
      if (Dune::Fem::Capabilities::isAdaptiveIndexSet<typename GridPart::IndexSetType>::v)
      {
        typedef VirtualizedGridFunction<GridPart,Dune::FieldVector<double,1>> Indicator;
        cls.def("mark",[](GridView &self,
                       Indicator &indicator,
                       double refineTolerance, double coarsenTolerance,
                       int minLevel, int maxLevel) {
            typename GridView::Grid &grid = const_cast<typename GridView::Grid&>(self.grid());
            return Dune::Fem::GridAdaptation::
                      mark( grid, indicator, refineTolerance, coarsenTolerance, minLevel, maxLevel );
          });
        cls.def("markNeighbors",[](GridView &self,
                       Indicator &indicator,
                       double refineTolerance, double coarsenTolerance,
                       int minLevel, int maxLevel) {
            typename GridView::Grid &grid = const_cast<typename GridView::Grid&>(self.grid());
            return Dune::Fem::GridAdaptation::
                      mark( grid, indicator, refineTolerance, coarsenTolerance, minLevel, maxLevel, true );
          });
        cls.def("doerflerMark",[](GridView &self, Indicator &indicator, double tolerance, int maxLevel,
                                  double layered) {
            typename GridView::Grid &grid = const_cast<typename GridView::Grid&>(self.grid());
            if (layered > 0)
              return Dune::Fem::GridAdaptation::layeredDoerflerMarking( grid, indicator, tolerance, maxLevel, layered );
            else
              return Dune::Fem::GridAdaptation::doerflerMarking( grid, indicator, tolerance, maxLevel );
          });
      }
    }
  }
}
#endif // DUNE_FEMPY_PY_GRIDVIEW
