#ifndef DUNE_FEM_MARKING_MAXIMUM_HH
#define DUNE_FEM_MARKING_MAXIMUM_HH

#include <cstddef>

#include <algorithm>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/function/common/discretefunction.hh>

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
      template <class Grid, class Indicator>
      static inline
      std::pair<int, int>
      doerflerMarking( Grid& grid, const Indicator& indicator,
                       const double theta, int maxLevel = -1 )
      {
        using std::max;

        typedef Dune::ReferenceElements< typename Grid::ctype, Grid::dimension > ReferenceElements;
        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;

        double totalError( 0 ), maxError( 0 );
        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;
          localIndicator.bind(e);
          const auto &center = ReferenceElements::general( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = std::abs(value[0]);
          totalError += eta;
          maxError = max( maxError, eta );
        }
        maxError = grid.comm().max( maxError );
        totalError = grid.comm().sum( totalError );

        // Let a.first be a real number and denote by a.second the sum of all
        // local errors greater than a.first.
        // Now consider two such pairs, a and b, such that
        //     a.second >= theta*totalError > b.second.
        // We seek to minimize this interval using bisection.
        std::pair< double, double > a( 0, totalError ), b( maxError, 0 );
        double m = (a.first + b.first) / double( 2 );
        while( (m > a.first) && (a.second > theta*totalError) )
        {
          // c.first: maximum local error less or equal to m
          // c.second: sum of local errors greater than m
          std::pair< double, double > c( 0, 0 );
          for (const auto &e : indicator.space())
          {
            if (!e.isLeaf()) continue;
            localIndicator.bind(e);
            const auto &center = ReferenceElements::general( e.type() ).position(0,0);
            localIndicator.evaluate(center,value);
            double eta = value[0];
            if( eta > m )
              c.second += eta;
            else
              c.first = max( c.first, eta );
          }
          c.first = grid.comm().max( c.first );
          c.second = grid.comm().sum( m );

          if( a.second > c.second )
          {
            // note: There is no local error in (m.first, c]
            if( c.second < theta*totalError )
              b = c;
            else
              a = std::make_pair( m, c.second );
            m = (a.first + b.first) / double( 2 );
          }
          else
          {
            // The total error could not be reduced. Is it still possible?
            // Try to find a local error m in (a.first, b.first) and use it as
            // next interval splitting point.
            // Note: If such an m exists, it must be greater than the mid point,
            //       because we already tried that one.
            m = 0;
            for (const auto &e : indicator.space())
            {
              if (!e.isLeaf()) continue;
              localIndicator.bind(e);
              const auto &center = ReferenceElements::general( e.type() ).position(0,0);
              localIndicator.evaluate(center,value);
              double eta = value[0];
              if( eta < b.first )
                m = max( m, eta );
            }
            m = grid.comm().max( m );
          }
        }

        // marking all elements with local error > a.first now yields the desired
        // property.
        std::size_t marked = 0;
        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;

          localIndicator.bind(e);
          const auto &center = ReferenceElements::general( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = value[0];
          if( eta <= a.first )
            continue;
          const auto &gridEntity = Dune::Fem::gridEntity(e);
          if (e.level()<maxLevel || maxLevel==-1)
          {
            grid.mark(Marker::refine, gridEntity);
            ++marked;
          }
          else
            grid.mark(Marker::keep, gridEntity);
        }

        return std::make_pair(marked,0);
      }


      /** \brief Layered Doerfler marking
       *  http://www.math.umd.edu/~rhn/teaching/m714/matlab/lecture-4.pdf
       */

      /** TODO: should be implemented without use of iterators, i.e., using
       * the dofvector directly
       **/
      namespace detail
      {
        template <class Grid, class Indicator>
        static inline std::pair< double, double >
        computeSumAndMaxGridWalk( Grid& grid, const Indicator& indicator,
                                  const double nu, std::vector< double >& buckets )
        {
          typedef Dune::ReferenceElements< typename Grid::ctype, Grid::dimension > ReferenceElements;
          Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
          typename Indicator::RangeType value;
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

          // compute global values
          maxIndicator = indicator.space().gridPart().comm().max( maxIndicator );
          sumIndicator = indicator.space().gridPart().comm().sum( sumIndicator );

          for (const auto &e : indicator.space())
          {
            if (!e.isLeaf()) continue;
            localIndicator.bind(e);
            const auto &center = ReferenceElements::general( e.type() ).position(0,0);
            localIndicator.evaluate(center,value);
            double eta = value[0];
            int index = int(eta/maxIndicator*1./nu);
            // std::cout << "   " << eta << " " << maxIndicator << " " << index << std::endl;
            assert( index < buckets.size() );
            buckets[index] += eta;
          }

          // compute global sum of all buckets in parallel
          indicator.space().gridPart().comm().sum( buckets.data(), buckets.size() );

          return std::make_pair( sumIndicator, maxIndicator );
        }

        template <class Grid, class Indicator>
        static inline std::pair< double, double >
        computeSumAndMax( Grid& grid, const Indicator& indicator,
                          const double nu, std::vector< double >& buckets )
        {
          return computeSumAndMaxGridWalk( grid, indicator, nu, buckets );
        }

        template <class Grid, class Imp>
        static inline std::pair< double, double >
        computeSumAndMax( Grid& grid, const Dune::Fem::DiscreteFunctionInterface< Imp >& indicator,
                          const double nu, std::vector< double >& buckets )
        {
          // if space is for some reason not constant then default to general method
          if( indicator.space().order() > 0 )
            return computeSumAndMaxGridWalk( grid, indicator, nu, buckets );

          double maxIndicator = 0;
          double sumIndicator = 0;
          // but this way we can avoid grid traversal etc.
          // at the moment we get a VirtualizedGF so this wouldn't work
          for (const auto &d : Dune::Fem::dofs(indicator) )
          {
            maxIndicator = std::max(maxIndicator,d);
            sumIndicator += d;
          }

          // compute global values
          maxIndicator = indicator.space().gridPart().comm().max( maxIndicator );
          sumIndicator = indicator.space().gridPart().comm().sum( sumIndicator );

          // let's assume that indicator is a FV function (would be good to be able to check this)
          // but this way we can avoid grid traversal etc.
          // at the moment we get a VirtualizedGF so this wouldn't work
          for (const auto &d : Dune::Fem::dofs(indicator) )
          {
            int index = int(d/maxIndicator*1./nu);
            assert( index < buckets.size() );
            buckets[index] += d;
          }

          // compute global sum of all buckets in parallel
          indicator.space().gridPart().comm().sum( buckets.data(), buckets.size() );

          return std::make_pair( sumIndicator, maxIndicator );
        }
      }

      template <class Grid, class Indicator>
      static inline
      std::pair<int, int>
      layeredDoerflerMarking( Grid& grid, const Indicator& indicator,
                              const double tolerance, int maxLevel = -1,
                              double nu = 0.05)
      {
        // if maxLevel < 0 set to maximum value
        maxLevel = ( maxLevel < 0 ) ? std::numeric_limits< int >::max() : maxLevel;

        int refMarked = 0;
        int crsMarked = 0;
        typedef Dune::ReferenceElements< typename Grid::ctype, Grid::dimension > ReferenceElements;
        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;

        // Part 1: compute sum and maximum of indicators
        // Part 2: subdivide [0,maxEta] into nu sized intervals and compute how
        // much of the overal error is in each interval
        std::vector<double> buckets(std::ceil(1./nu)+1);
        auto sumMax = detail::computeSumAndMax( grid, indicator, nu, buckets );
        double sumIndicator = sumMax.first;
        double maxIndicator = sumMax.second;

        // Part 3: compute how many buckets we have to use to get
        // above (1-tolerance)*(1-tolerance)*sumIndicator:
        double gamma = 1;
        double sum = 0;
        for (int index=buckets.size()-1;
             index >= 0 && sum < (1-tolerance)*(1-tolerance)*sumIndicator;
             --index)
        {
          sum += buckets[index];
          gamma -= nu;
        }

        // no communication should be required and we can simply now
        // go ahead with the refinement:
        //
        sum = 0; // just for checkint that it worked...
        const double gammaMaxIndicator = gamma*maxIndicator ;
        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;
          const auto &gridEntity = Dune::Fem::gridEntity(e);
          localIndicator.bind(e);
          const auto &center = ReferenceElements::general( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = value[0];
          if (eta > gammaMaxIndicator )
          {
            if (e.level()<maxLevel)
              refMarked += grid.mark(Marker::refine, gridEntity);
            else
              grid.mark(Marker::keep, gridEntity);
            // although we might not have marked for refinement due to
            // level restriction we count this as part of the indicator
            // taken care of
            sum += eta;
          }
        }
        // just checking that the algorihtm did the right thing...
        assert( sum >= (1-tolerance)*(1-tolerance)*sumIndicator);
        return std::make_pair( grid.comm().sum(refMarked), grid.comm().sum(crsMarked) );
      }

    } // end namespace GridAdaptation

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_MAXIMUM_HH
