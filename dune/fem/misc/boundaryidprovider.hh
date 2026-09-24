#ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
#define DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typeutilities.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/function/common/localcontribution.hh>

namespace Dune
{
  namespace Fem
  {

    template< class Grid >
    struct BoundaryIdProvider
    {
      template <class T, class = int>
      struct hasBoundaryId : std::false_type {};

      template <class T>
      struct hasBoundaryId<T, decltype(std::declval<T>().impl().boundaryId())> : std::true_type {};

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        static constexpr bool hasBndId = hasBoundaryId< Intersection >::value;

        // for all grids that have the method, call it on the implementation
        if constexpr ( hasBndId )
        {
          return intersection.impl().boundaryId();
        }
        else // otherwise use fallback
        {
          // Cartesian grids use indexInInside + 1
          if constexpr ( Dune::Capabilities::isCartesian< Grid > :: v )
          {
            return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
          }
          else // e.g. UGGrid
          {
            return intersection.boundarySegmentIndex();
          }
        }
      }
    };

    // BoundaryIdProvider for general GridParts or GridViews
    // -----------------------------------------------------

    //! this works for both, GridView and GridPart
    template< class GridView, class Intersection>
    inline static int boundaryId ( const Intersection &intersection )
    {
      return Dune::Fem::BoundaryIdProvider< typename GridView::Grid > ::
             boundaryId( intersection );
    }

    //! this works for both, GridView and GridPart
    template< class GridView, class Intersection>
    inline static int boundaryId ( const GridView&, const Intersection &intersection )
    {
      return Dune::Fem::BoundaryIdProvider< typename GridView::Grid > ::
             boundaryId( intersection );
    }


    /* \brief BoundaryId inspection: project boundary ids to piecewise constant
              function. Elements with more than one boundary segment will
              contain the sum of all boundary ids.
    */
    template <class DiscreteFunction>
    void projectBoundaryIds( DiscreteFunction& df )
    {
      if( df.space().order() != 0 ) // should be piecewise constant, and FV space
      {
        DUNE_THROW(InvalidStateException,"projectBoundaryIds: expect piecewise constant discrete function");
      }

      df.clear(); // reset all dofs

      const auto& gridPart = df.space().gridPart();

      typedef AddLocalContribution< DiscreteFunction > AddLocalContributionType;
      AddLocalContributionType localDf( df );
      for( const auto& element : df.space() )
      {
        int count = 0;
        auto guard = bindGuard( localDf, element );
        for( const auto& intersection : intersections( gridPart, element ) )
        {
          if( intersection.boundary() )
          {
            localDf[ 0 ] += boundaryId( gridPart, intersection );
            count += 1;
          }
        }
        if (count>0)
          localDf[ 0 ] /= double(count);
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
