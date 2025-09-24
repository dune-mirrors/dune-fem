#ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
#define DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
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
      typedef typename GridPartType::GridType GridType;
      typedef DofManager< GridType > DofManagerType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename EntityType::Geometry::GlobalCoordinate        GlobalCoordinateType;

      //typedef FunctionSpace< double, int, GridPartType::dimensionworld, 1 >  FunctionSpaceType;
      // typedef FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0, SimpleStorage > SpaceType;

      //typedef ParallelDofMapper< GridPartType, typename SpaceType::BlockMapperType > ParallelMapperType;
      //typedef typename ParallelMapperType :: GlobalKeyType  GlobalKeyType;

      static const int dimensionworld = GridPartType::dimensionworld;

      explicit DefaultUniqueFacetOrientation ( const GridPartType &gridPart )
        : gridPart_( gridPart ),
          dofManager_( DofManagerType::instance( gridPart_.grid() ) ),
          uniqueDirection_( M_PI ),
          sequence_( -1 )
      {
        if constexpr ( dimensionworld > 1 )
          uniqueDirection_[ 1 ] = M_LN2;
        if constexpr ( dimensionworld > 2 )
          uniqueDirection_[ 2 ] = M_E;
      }

      void update() const
      {
        if( sequence_ != dofManager_.sequence() )
        {
          orientations_.resize( gridPart().indexSet().size(0) );
          const auto& indexSet = gridPart().indexSet();
          for( const auto& entity : elements( gridPart(), Partitions::all ) )
          {
            unsigned int orientations = 0;
            for( auto intersection : intersections( gridPart(), entity ) )
            {
              if( intersection.neighbor() )
              {
                // make sure that the normal is not orthogonal to the uniqueDirection_
                assert( std::abs(intersection.centerUnitOuterNormal() * uniqueDirection_) > 0 );

                if( (intersection.centerUnitOuterNormal() * uniqueDirection_) > 0 )
                {
                  orientations |= (1 << intersection.indexInInside());
                }
              }
            }
            orientations_[ indexSet.index( entity ) ] = orientations;
          }
          sequence_ = dofManager_.sequence();
        }
      }

      unsigned int operator() ( const EntityType &entity ) const
      {
        update();
        return orientations_[ gridPart().indexSet().index( entity ) ];
      }

        /*
        unsigned int orientations = 0;
        for( auto intersection : intersections( gridPart(), entity ) )
        {
          if( intersection.neighbor() && (globallyUniqueIndex( entity ) < globallyUniqueIndex( intersection.outside() )) )
            orientations |= (1 << intersection.indexInInside());
        }
        return orientations;
        */

      /*
      template <class Entity>
      GlobalKeyType globallyUniqueIndex( const Entity& entity ) const
      {
        retu
        GlobalKeyType index = -1;
        // parallelMapper_.mapEach( entity, [ &index ] ( auto local, auto global ) { assert( local == 0 ); index = global; } );
        return gridPart().indexSet().index(entity);
        return index;
      }
      */

      const GridPartType &gridPart () const { return gridPart_; }

    protected:
      const GridPartType& gridPart_;
      const DofManagerType& dofManager_;
      GlobalCoordinateType uniqueDirection_;
      mutable std::vector< unsigned int > orientations_;
      //SpaceType           space_;  // works without this
      // mutable ParallelMapperType  parallelMapper_;
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

    template< class GP >
    using UniqueFacetOrientation =
      typename std::conditional< GridPartCapabilities::isCartesian< GP >::v,
                                 CartesianUniqueFacetOrientation< GP>,
                                 DefaultUniqueFacetOrientation< GP > >::type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_UNIQUEFACETORIENTATION_HH
