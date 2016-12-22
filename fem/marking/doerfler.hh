#ifndef DUNE_FEM_MARKING_MAXIMUM_HH
#define DUNE_FEM_MARKING_MAXIMUM_HH

#include <cstddef>

#include <algorithm>

#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  namespace Fem
  {

    // doerflerMarking
    // ---------------

    /** \brief doerflerMarking
     *
     * Mark a minimal set \f$\mathcal{A} \subset \mathcal{G}\f$ of elements in
     * a grid \f$\mathcal{G}\f$ such that
     * \$[
     *   \sum_{T \in \mathcal{A}} \eta_t \ge \theta\,\sum_{T \in \mathcal{G}} \eta_T.
     * \$]
     *
     * See also:
     * W. Dörfler, A Convergent Adaptive Algorithm for Poisson's Equation,
     * SIAM J. Numer. Anal. 33 (3), 1106-1124, 1996
     *
     * For the sake of simplicity, this algorithm assumes disjoint local errors
     * \f$\eta_T\f$. Otherwise, too more elements may be marked.
     *
     * \param[in]     localError  function modelling \f$\eta_T\f$
     * \param[in]     theta       factor of total error to mark
     * \param[inout]  grid        grid \f$\mathcal{G}\f$ to mark
     **/
    template< class LocalError, class Real, class Grid >
    inline static std::size_t doerflerMarking ( const LocalError &localError, const Real &theta, Grid &grid )
    {
      using std::max;

      const auto gridView = grid.leafGridView();

      Real totalError( 0 ), maxError( 0 );
      for( const auto &element : elements( grid.leafGridView(), Partitions::interior ) )
      {
        const Real e = localError( element );
        totalError += e;
        maxError = max( maxError, e );
      }
      maxError = grid.comm().max( maxError );
      totalError = gridView.comm().sum( totalError );

      // Let a.first be a real number and denote by a.second the sum of all
      // local errors greater than a.first.
      // Now consider two such pairs, a and b, such that
      //     a.second >= theta*totalError > b.second.
      // We seek to minimize this interval using bisection.
      std::pair< Real, Real > a( 0, totalError ), b( maxError, 0 );
      Real m = (a.first + b.first) / Real( 2 );
      while( (m > a.first) && (a.second > theta*totalError) )
      {
        // c.first: maximum local error less or equal to m
        // c.second: sum of local errors greater than m
        std::pair< Real, Real > c( 0, 0 );
        for( const auto &element : elements( gridView, Partitions::interior ) )
        {
          const Real e = localError( element );
          if( e > m )
            c.second += e;
          else
            c.first = max( c.first, e );
        }
        c.first = gridView.comm().max( c.first );
        c.second = gridView.comm().sum( m.second );

        if( a.second > c.second )
        {
          // note: There is no local error in (m.first, c]
          if( c.second < theta*totalError )
            b = c;
          else
            a = std::make_pair( m, c.second );
          m = (a.first + b.first) / Real( 2 );
        }
        else
        {
          // The total error could not be reduced. Is it still possible?
          // Try to find a local error m in (a.first, b.first) and use it as
          // next interval splitting point.
          // Note: If such an m exists, it must be greater than the mid point,
          //       because we already tried that one.
          m = 0;
          for( const auto &element : elements( gridView, Partitions::interior ) )
          {
            const Real e = localError( element );
            if( e < b.first )
              m = max( m, e );
          }
          m = gridView.comm.max( m );
        }
      }

      // marking all elements with local error > a.first now yields the desired
      // property.
      std::size_t marked = 0;
      for( const auto &element : elements( grid.leafGridView() ) )
      {
        if( localError( element ) <= a.first )
          continue;
        grid_.mark( 1, element );
        ++marked;
      }
      return grid.comm().sum( marked );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
