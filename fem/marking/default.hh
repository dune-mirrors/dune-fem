#ifndef DUNE_FEMPY_PY_DEFAULT_MARKING_HH
#define DUNE_FEMPY_PY_DEFAULT_MARKING_HH

#include <memory>
#include <dune/geometry/referenceelements.hh>
#include <dune/fempy/py/grid/gridpart.hh>

namespace Dune
{
  namespace Fem
  {
    namespace GridAdaptation
    {
      struct Marker
      {
        static const int refine  =  1;
        static const int keep    =  0;
        static const int coarsen = -1;
      };

      // Indicator is of type DiscreteFunction
      template <class Grid, class Indicator>
      static inline
      std::pair< int, int >
      mark( Grid& grid, Indicator& indicator,
            const double refineTolerance, const double coarsenTolerance,
            int minLevel = 0, int maxLevel = -1 )
      {
        typedef Dune::ReferenceElements< typename Grid::ctype, Grid::dimension > ReferenceElements;

        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;

        int refMarked = 0;
        int crsMarked = 0;

        if ( maxLevel < 0 )
          maxLevel = std::numeric_limits<int>::max();

        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;

          const auto &gridEntity = Dune::Fem::gridEntity(e);
          int marked = grid.getMark(gridEntity);
          if (marked==1) continue;

          localIndicator.bind(e);
          const auto &center = ReferenceElements::general( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = value[0];
          const int level = e.level();
          if (eta>refineTolerance && level<maxLevel)
            refMarked += grid.mark( Marker::refine, gridEntity);
          else if (eta<coarsenTolerance && level>minLevel)
            crsMarked += grid.mark( Marker::coarsen, gridEntity);
          else
            grid.mark(Marker::keep, gridEntity);
        }
        return std::make_pair( grid.comm().sum(refMarked), grid.comm().sum(crsMarked) );
      } // end mark

    } // end namespace GridAdaptation

  } // end namespace Fem
} // end namespace Dune
#endif
