#ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
#define DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/mapper/parallel.hh>

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

      typedef FunctionSpace< double, int, GridPartType::dimensionworld, 1 >  FunctionSpaceType;
      typedef FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0, SimpleStorage > SpaceType;

      typedef ParallelDofMapper< GridPartType, typename SpaceType::BlockMapperType > ParallelMapperType;
      typedef typename ParallelMapperType :: GlobalKeyType  GlobalKeyType;

      explicit DefaultUniqueFacetOrientation ( const GridPartType &gridPart )
        : gridPart_( gridPart ),
          space_( const_cast< GridPartType& > (gridPart_) ),
          parallelMapper_( gridPart_, space_.blockMapper() ),
          sequence_( -1 )
      {
      }

      unsigned int operator() ( const EntityType &entity ) const
      {
        if( sequence_ != space_.sequence() )
        {
          parallelMapper_.update();
          sequence_ = space_.sequence();
        }

        unsigned int orientations = 0;
        for( auto intersection : intersections( gridPart(), entity ) )
        {
          if( intersection.neighbor() && (globallyUniqueIndex( entity ) < globallyUniqueIndex( intersection.outside() )) )
            orientations |= (1 << intersection.indexInInside());
        }
        return orientations;
      }

      template <class Entity>
      GlobalKeyType globallyUniqueIndex( const Entity& entity ) const
      {
        GlobalKeyType index = -1;
        parallelMapper_.mapEach( entity, [ &index ] ( auto local, auto global ) { assert( local == 0 ); index = global; } );
        return index;
      }

      const GridPartType &gridPart () const { return gridPart_; }

    protected:
      const GridPartType& gridPart_;
      SpaceType           space_;
      mutable ParallelMapperType  parallelMapper_;
      mutable int sequence_;
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

    template< class GridPart >
    using UniqueFacetOrientation =
      typename std::conditional< GridPartCapabilities::isCartesian< GridPart >::v,
                                 CartesianUniqueFacetOrientation< GridPart>,
                                 DefaultUniqueFacetOrientation< GridPart > >::type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
