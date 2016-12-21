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

    template< class LocalError, class Real, class Grid >
    inline static std::size_t maximalMarking ( const LocalError &localError, const Real &tolerance, Grid &grid )
    {
      using std::max;

      Real maxError( 0 );
      for( const auto &element : elements( grid.leafGridView() ) )
        maxError = max( maxError, localError( element ) );
      maxError = grid.comm().max( maxError );

      std::size_t marked = 0;
      for( const auto &element : elements( grid.leafGridView() ) )
      {
        if( localError( element ) <= tolerance*maxError )
          continue;
        grid_.mark( 1, element );
        ++marked;
      }
      return grid.comm().sum( marked );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
