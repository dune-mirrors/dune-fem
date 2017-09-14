#ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
#define DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // DefaultUniqueFacetOrientation
    // -----------------------------

    template< class GridPart >
    struct DefaultUniqueFacetOrientation
    {
      typedef GridPart GridPartType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      explicit DefaultUniqueFacetOrientation ( const GridPartType &gridPart )
        : gridPart_( gridPart )
      {}

      unsigned int operator() ( const EntityType &entity ) const
      {
        unsigned int orientations = 0;
        const auto &indexSet = gridPart().indexSet();
        for( auto intersection : intersections( gridPart(), entity ) )
        {
          if( intersection.neighbor() && (indexSet.index( entity ) < indexSet.index( intersection.outside() )) )
            orientations |= (1 << intersection.indexInInside());
        }
        return orientations;
      }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      const GridPartType &gridPart_;
    };



    // CartesianUniqueFacetOrientation
    // -------------------------------

    template< class GridPart >
    struct CartesianUniqueFacetOrientation
    {
      typedef GridPart GridPartType;
      typedef typename GridPart::template Codim< 0 >::EntityType EntityType;

      explicit CartesianUniqueFacetOrientation ( const GridPartType &gridPart )
        : gridPart_( gridPart )
      {}

      unsigned int operator() ( const EntityType &entity ) const
      {
        return orientations( std::make_index_sequence< GridPartType::dimension >() );
      }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      template< std::size_t... i >
      static constexpr unsigned int orientations ( std::index_sequence< i... > )
      {
        return Std::sum( (1u << 2*i)... );
      }

      const GridPartType &gridPart_;
    };



    // UniqueFacetOrientation
    // ----------------------

    namespace Impl
    {

      template< class GridPart, class = void >
      struct UniqueFacetOrientationType
      {
        typedef DefaultUniqueFacetOrientation< GridPart > Type;
      };

      template< class GridPart >
      struct UniqueFacetOrientationType< GridPart, std::enable_if_t< GridPartCapabilities::isCartesian< GridPart >::v > >
      {
        typedef CartesianUniqueFacetOrientation< GridPart > Type;
      };

    } // namespace Impl

    template< class GridPart >
    using UniqueFacetOrientation = typename Impl::UniqueFacetOrientationType< GridPart >::Type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
