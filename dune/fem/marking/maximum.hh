#ifndef DUNE_FEM_MARKING_MAXIMUM_HH
#define DUNE_FEM_MARKING_MAXIMUM_HH

#include <cstddef>

#include <algorithm>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/function/localfunction/const.hh>

#include <dune/fem/marking/default.hh>

namespace Dune
{

  namespace Fem
  {

    // maximalMarking
    // --------------

    template <class Grid, class Indicator>
    static inline
    std::pair<int, int>
    maximumMarking( Grid& grid, const Indicator& indicator,
                    double theta, int maxLevel = -1)
    {
      using Marker = GridAdaptation::Marker;

      int refMarked = 0;
      int crsMarked = 0;
      maxLevel = ( maxLevel < 0 ) ? std::numeric_limits< int >::max() : maxLevel;
      using std::max;

      double maxError( 0 );
      for( const auto &element : elements( grid.leafGridView() ) )
        maxError = max( maxError, localError( element ) );
      maxError = grid.comm().max( maxError );

      ConstLocalFunction< Indicator > localIndicator( indicator );
      typename Indicator :: RangeType value;

      for (const auto &e : indicator.space())
      {
        if (!e.isLeaf()) continue;
        const auto &gridEntity = Dune::Fem::gridEntity(e);
        localIndicator.bind(e);
        const auto& center = Dune::referenceElement< typename Grid::ctype, Grid::dimension>( e.type() ).position(0,0);
        localIndicator.evaluate(center,value);
        double eta = value[0];
        if( eta > theta*maxError )
        {
          if (e.level()<maxLevel)
            refMarked += grid.mark(Marker::refine, gridEntity);
          else
            grid.mark(Marker::keep, gridEntity);
        }
      }
      return std::make_pair( grid.comm().sum(refMarked), grid.comm().sum(crsMarked) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
