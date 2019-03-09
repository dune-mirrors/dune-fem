#ifndef DUNE_FEM_MARKING_MAXIMUM_HH
#define DUNE_FEM_MARKING_MAXIMUM_HH

#include <cstddef>

#include <algorithm>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/marking/default.hh>

namespace Dune
{

  namespace Fem
  {

    namespace GridAdaptation
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
        for( const auto &element : elements( gridView, Partitions::interior ) )
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
        for( const auto &element : elements( gridView, Partitions::interior ) )
        {
          if( localError( element ) <= a.first )
            continue;
          grid.mark( Marker::refine, element );
          ++marked;
        }
        return grid.comm().sum( marked );
      }


      /** \brief Layered Doerfler marking
       *  http://www.math.umd.edu/~rhn/teaching/m714/matlab/lecture-4.pdf
       */

      template <class Grid, class Indicator>
      static inline
      std::pair<int, int>
      layeredDoerflerMarking( Grid& grid, const Indicator& indicator,
                              const double tolerance, int maxLevel = -1 )
      {
        typedef Dune::ReferenceElements< typename Grid::ctype, Grid::dimension > ReferenceElements;
        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;
        int refMarked = 0;
        double maxIndicator = 0;
        double sumIndicator = 0;
        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;
          localIndicator.bind(e);
          const auto &center = ReferenceElements::general( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = value[0];
          maxIndicator = std::max(maxIndicator,eta);
          sumIndicator += eta;
        }
        double gamma = 1;
        double nu = 0.05;
        double sum = 0;
        double oldSum = -1;
        while ( sum < (1-tolerance)*(1-tolerance)*sumIndicator )
        {
          gamma -= nu;
          for (const auto &e : indicator.space())
          {
            if (!e.isLeaf()) continue;

            const auto &gridEntity = Dune::Fem::gridEntity(e);
            int marked = grid.getMark( gridEntity );
            // if (marked==1) continue;

            localIndicator.bind(e);
            const auto &center = ReferenceElements::general( e.type() ).position(0,0);
            localIndicator.evaluate(center,value);
            double eta = value[0];
            if (eta > gamma*sumIndicator && (e.level()<maxLevel || maxLevel==-1))
            {
              refMarked += grid.mark(Marker::refine, gridEntity);
              sum += eta;
            }
            else
              grid.mark(Marker::keep, gridEntity);
          }
          if (std::abs(oldSum-sum)<1e-15) break;
        }
        return std::make_pair(refMarked,0);
      }

    } // end namespace GridAdaptation

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
