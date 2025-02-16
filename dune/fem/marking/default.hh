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

      // Indicator is of type DiscreteFunction or GridFunction
      template <class Grid, class Indicator>
      static inline
      std::pair< int, int >
      mark( Grid& grid, Indicator& indicator,
            const double refineTolerance, const double coarsenTolerance,
            const int minLevel = 0,
            int maxLevel = -1,
            const double minVolume = -1.,
            double maxVolume = -1.0,
            const bool markNeighbors = false,
            const bool statistics = false ) // if true return number of marked cells
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
          const double volume = (useVolumeRefinement || useVolumeCoarsening) ? localIndicator.geometry().volume() : 0.0;

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

        // only do communication if statistics is enabled
        if( statistics )
        {
          gridPart.comm().sum( sumMarks.data(), 2 );
          return std::make_pair( refMarked, crsMarked );
        }
        else
          return std::make_pair( int(-1), int(-1));
      } // end mark

    } // end namespace GridAdaptation


    namespace SpaceAdaptation
    {
      // Indicator is of type DiscreteFunction or GridFunction
      template <class Space, class Indicator>
      static inline
      std::pair< int, int >
      mark( Space& space, Indicator& indicator,
            const double refineTolerance, const double coarsenTolerance, // tolerances for increasing or decreasing polynomial order
            const int minOrder, const int maxOrder, // min and max polOrders [0, space.order()]
            const bool markNeighbors = false,  // also mark neighbors
            const bool statistics = false ) // if true return number of marked cells
      {
        Dune::Fem::ConstLocalFunction<Indicator> localIndicator(indicator);
        typename Indicator::RangeType value;
        if( minOrder < 0 && minOrder > space.order() )
          DUNE_THROW(Dune::InvalidStateException,"SpaceAdaptation::mark: minOrder out of range: select something between 0 and the space's order! Selected was " << minOrder);

        if( maxOrder < minOrder && maxOrder > space.order() )
          DUNE_THROW(Dune::InvalidStateException,"SpaceAdaptation::mark: maxOrder out of range: select something between minOrder and the space's order! Selected was " << maxOrder);

        const auto& gridPart = space.gridPart();

        std::array< int, 2 > sumMarks = {{ 0,0 }};
        int& refMarked = sumMarks[ 0 ]; // increased polynomial order
        int& crsMarked = sumMarks[ 1 ]; // decreased polynomial order

        std::vector< bool > marked;
        if( markNeighbors )
          marked.resize( gridPart.indexSet().size( 0 ), false );

        for (const auto &e : space)
        {
          localIndicator.bind(e);
          const auto &center = Dune::referenceElement( e ).position(0,0);
          localIndicator.evaluate(center,value);
          double eta = std::abs(value[0]);

          const bool markedAlready = ( markNeighbors ) ? marked[ gridPart.indexSet().index( e ) ] : false;

          // get current polynomial degree
          const int currentOrder = space.order( e );
          int polOrder = currentOrder;

          // if the element itself has been marked check that
          // the current mark is lower than the proposed mark
          if( markedAlready )
          {
            // if the current mark is higher than the current order just continue
            if( space.getMark( e ) > polOrder )
              continue ;
          }

          // check that estimator is larger than tolerance
          // check that polOrder is larger than minimal polOrder
          if (eta > refineTolerance && polOrder < maxOrder)
          {
            // increase order by one
            ++ polOrder;

            if( markNeighbors )
            {
              const auto iEnd = gridPart.iend( e );
              for( auto it = gridPart.ibegin( e ); it != iEnd; ++it )
              {
                const auto& intersection = *it;
                if( intersection.neighbor() )
                {
                  const auto& outside = intersection.outside();
                  // mark neighbor for refinement only (i.e. if polynomial order is higher)
                  const int nbOrder = space.getMark( outside );
                  if ( polOrder > nbOrder )
                  {
                    space.mark( polOrder, outside );
                    ++refMarked;

                    // mark as marked (marked should be resized correctly)
                    marked[ gridPart.indexSet().index( outside ) ] = true;
                  }
                }
              }
            }
          }
          else if (! markedAlready && eta < coarsenTolerance && polOrder > minOrder)
          {
            -- polOrder;
          }
          // else keep the same polynomial order

          // increase counters
          if( polOrder > currentOrder )
            ++refMarked;
          else if ( polOrder < currentOrder )
            ++crsMarked;

          // mark polynomial order for current element
          space.mark( polOrder, e );

          // mark as marked
          if( markNeighbors )
            marked[ gridPart.indexSet().index( e ) ] = true;
        }

        // only do communication if statistics is enabled
        if( statistics )
        {
          gridPart.comm().sum( sumMarks.data(), 2 );
          return std::make_pair( refMarked, crsMarked );
        }
        else
          return std::make_pair( int(-1), int(-1));
      } // end mark

    } // end namespace SpaceAdaptation

  } // end namespace Fem
} // end namespace Dune
#endif
