#ifndef DUNE_FEMPY_PY_GRIDVIEW
#define DUNE_FEMPY_PY_GRIDVIEW
#include <dune/python/grid/gridview.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>
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
            Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
            typename Indicator::RangeType value;
            int refMarked = 0;
            int crsMarked = 0;
            if (maxLevel==-1)
              maxLevel = std::numeric_limits<int>::max();
            for (const auto &e : indicator.space())
            {
              if (!e.isLeaf()) continue;
              const auto &entity = Dune::Fem::gridEntity(e);
              localIndicator.bind(e);
              const auto &center = ReferenceElements< typename GridView::ctype, GridView::dimension >::general( e.type() ).position(0,0);
              localIndicator.evaluate(center,value);
              double eta = value[0];
              int marked = grid.getMark(entity);
              if (marked==1) continue;
              if (eta>refineTolerance && e.level()<maxLevel)
                refMarked += grid.mark(1, entity);
              else if (eta<coarsenTolerance && e.level()>minLevel)
                crsMarked += grid.mark(-1, entity);
              else
                grid.mark(0, entity);
            }
            return std::make_pair(refMarked,crsMarked);
          });
        cls.def("doerflerMark",[](GridView &self,
                       Indicator &indicator,
                       double tolerance, int maxLevel) {
            // layered doerfler strategy
            // http://www.math.umd.edu/~rhn/teaching/m714/matlab/lecture-4.pdf
            if (maxLevel==-1)
              maxLevel = std::numeric_limits<int>::max();
            typename GridView::Grid &grid = const_cast<typename GridView::Grid&>(self.grid());
            Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
            typename Indicator::RangeType value;
            int refMarked = 0;
            double maxIndicator = 0;
            double sumIndicator = 0;
            for (const auto &e : indicator.space())
            {
              if (!e.isLeaf()) continue;
              const auto &entity = Dune::Fem::gridEntity(e);
              localIndicator.bind(e);
              const auto &center = ReferenceElements< typename GridView::ctype, GridView::dimension >::general( e.type() ).position(0,0);
              localIndicator.evaluate(center,value);
              double eta = value[0];
              maxIndicator = std::max(maxIndicator,eta);
              sumIndicator += eta;
            }
            double gamma = 1;
            double nu = 0.05;
            double sum = 0;
            while ( sum < (1-tolerance)*(1-tolerance)*sumIndicator )
            {
              gamma -= nu;
              for (const auto &e : indicator.space())
              {
                if (!e.isLeaf()) continue;
                const auto &entity = Dune::Fem::gridEntity(e);
                localIndicator.bind(e);
                const auto &center = ReferenceElements< typename GridView::ctype, GridView::dimension >::general( e.type() ).position(0,0);
                localIndicator.evaluate(center,value);
                double eta = value[0];
                int marked = grid.getMark(entity);
                if (marked==1) continue;
                if (eta > gamma*sumIndicator && e.level()<maxLevel)
                {
                  refMarked += grid.mark(1, entity);
                  sum += eta;
                }
                else
                  grid.mark(0, entity);
              }
            }
            return std::make_pair(refMarked,0);
          });
      }
    }
  }
}
#endif // DUNE_FEMPY_PY_GRIDVIEW
