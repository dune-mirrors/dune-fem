#ifndef DUNE_FEMPY_PY_DEFAULT_MARKING_HH
#define DUNE_FEMPY_PY_DEFAULT_MARKING_HH

#include <array>
#include <memory>
#include <dune/geometry/referenceelements.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/function/localfunction/const.hh>

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
            const int minLevel = 0,
            int maxLevel = -1,
            const double minVolume = -1.,
            double maxVolume = -1.0,
            const bool markNeighbors = false )
      {
        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;

        std::array< int, 2 > sumMarks = {{ 0,0 }};

        int& refMarked = sumMarks[ 0 ];
        int& crsMarked = sumMarks[ 1 ];

        if ( maxLevel < 0 )
          maxLevel = std::numeric_limits<int>::max();

        const bool useVolumeRefinement = minVolume > 0.0;
        const bool useVolumeCoarsening = maxVolume > 0.0;
        if ( ! useVolumeCoarsening )
        {
          // this will avoid volume check
          maxVolume = std::numeric_limits<double>::max();
        }

        const auto& gridPart = indicator.space().gridPart();

        for (const auto &e : indicator.space())
        {
          if (!e.isLeaf()) continue;

          const auto &gridEntity = Dune::Fem::gridEntity(e);
          int marked = grid.getMark(gridEntity);
          if (marked==1) continue;

          localIndicator.bind(e);
          const auto &center = Dune::referenceElement< typename Grid::ctype, Grid::dimension>( e.type() ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = std::abs(value[0]);
          const int level = e.level();

          // compute volume only if necessary
          const double volume = (useVolumeRefinement || useVolumeCoarsening) ? e.geometry().volume() : 0.0;

          // check that estimator is larger than tolerance
          // check that level is smaller than minimal level
          // check that volume of element is larger than minimal acceptable volume
          if (eta>refineTolerance && level<maxLevel && volume>minVolume)
          {
            refMarked += grid.mark( Marker::refine, gridEntity);
            if( markNeighbors )
            {
              const auto iEnd = gridPart.iend( e );
              for( auto it = gridPart.ibegin( e ); it != iEnd; ++it )
              {
                const auto& intersection = *it;
                if( intersection.neighbor() )
                {
                  const auto& outside = intersection.outside();
                  const double outsideVol = (useVolumeRefinement) ? outside.geometry().volume() : 0.0;
                  if (outside.level()<maxLevel && outsideVol>minVolume )
                    refMarked += grid.mark( Marker::refine, Dune::Fem::gridEntity(outside) );
                }
              }

            }
          }
          else if (eta<coarsenTolerance && level>minLevel && volume<maxVolume )
            crsMarked += grid.mark( Marker::coarsen, gridEntity);
          else
            grid.mark(Marker::keep, gridEntity);
        }

        grid.comm().sum( sumMarks.data(), 2 );
        return std::make_pair( refMarked, crsMarked );
      } // end mark

    } // end namespace GridAdaptation

  } // end namespace Fem
} // end namespace Dune
#endif
