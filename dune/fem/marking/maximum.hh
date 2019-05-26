#ifndef DUNE_FEM_MARKING_MAXIMUM_HH
#define DUNE_FEM_MARKING_MAXIMUM_HH

#include <cstddef>

#include <algorithm>

#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  namespace Fem
  {

    // maximalMarking
    // --------------

    static inline
    std::pair<int, int>
    maximumMarking( Grid& grid, const Indicator& indicator,
                    double theta, int maxLevel = -1)
    {
      int refMarked = 0;
      int crsMarked = 0;
      maxLevel = ( maxLevel < 0 ) ? std::numeric_limits< int >::max() : maxLevel;
      using std::max;

      double maxError( 0 );
      for( const auto &element : elements( grid.leafGridView() ) )
        maxError = max( maxError, localError( element ) );
      maxError = grid.comm().max( maxError );

      for (const auto &e : indicator.space())
      {
        if (!e.isLeaf()) continue;
        const auto &gridEntity = Dune::Fem::gridEntity(e);
        localIndicator.bind(e);
        const auto &center = ReferenceElements::general( e.type() ).position(0,0);
        localIndicator.evaluate(center,value);
        double eta = value[0];
        if( eta > theta*maxError )
          if (e.level()<maxLevel)
            refMarked += grid.mark(Marker::refine, gridEntity);
          else
            grid.mark(Marker::keep, gridEntity);
      }
      return std::make_pair( grid.comm().sum(refMarked), grid.comm().sum(crsMarked) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
