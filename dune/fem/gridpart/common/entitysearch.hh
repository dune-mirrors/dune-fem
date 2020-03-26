#ifndef DUNE_FEM_GRIDPART_ENTITYSEARCH_HH
#define DUNE_FEM_GRIDPART_ENTITYSEARCH_HH

#include <type_traits>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // DefaultEntitySearch
    // -------------------

    template< class GridPart, int codim, PartitionIteratorType partition >
    class DefaultEntitySearch
    {
      typedef DefaultEntitySearch< GridPart, codim, partition > ThisType;

      static const int dimension = GridPart::dimension;
      static const int dimensionworld = GridPart::dimensionworld;
      static const int codimension = codim;
      static const int mydimension = dimension - codimension;

      typedef typename GridPart::template Codim< codimension >::GeometryType GeometryType;
      typedef typename GridPart::template Codim< codimension >::template Partition< partition >::IteratorType IteratorType;

      typedef typename GeometryType::ctype ctype;
      typedef typename GeometryType::LocalCoordinate LocalCoordinateType;

    public:
      typedef GridPart GridPartType;

      typedef typename GridPart::template Codim< codimension >::EntityType EntityType;

      typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;

      explicit DefaultEntitySearch ( const GridPartType &gridPart )
      : gridPart_( gridPart )
      {}

      EntityType operator() ( const GlobalCoordinateType &x ) const
      {
        const auto end = gridPart_.template end< codimension, partition >();
        for( auto it = gridPart_.template begin< codimension, partition >(); it != end; ++it )
        {
          const auto& entity = *it;
          const auto geo = entity.geometry();
          const auto z = geo.local( x );
          if( (mydimension < dimensionworld) && ((geo.global( z ) - x).two_norm() > 1e-8 ) )
            continue;

          if( referenceElement<ctype,mydimension>( geo.type() ).checkInside( z ) )
            return entity;
        }
        DUNE_THROW( GridError, "Coordinate " << x << " is outside the grid." );
      }

    private:
      const GridPartType &gridPart_;
    };



    // GridEntitySearch
    // ----------------

    template< class GridPart, int codim, PartitionIteratorType partition >
    class GridEntitySearch
    : public DefaultEntitySearch< GridPart, codim, partition >
    {
      typedef GridEntitySearch< GridPart, codim, partition > ThisType;
      typedef DefaultEntitySearch< GridPart, codim, partition > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit GridEntitySearch ( const GridPartType &gridPart )
      : BaseType( gridPart )
      {}
    };

    template< class GridPart, PartitionIteratorType partition >
    class GridEntitySearch< GridPart, 0, partition >
    {
      typedef GridEntitySearch< GridPart, 0, partition > ThisType;

      static const int dimension = GridPart::dimension;
      static const int dimensionworld = GridPart::dimensionworld;
      static const int codimension = 0;
      static const int mydimension = dimension - codimension;

      typedef typename GridPart::template Codim< codimension >::GeometryType GeometryType;

    public:
      typedef GridPart GridPartType;

      typedef typename GridPart::template Codim< codimension >::EntityType EntityType;

      typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;

      explicit GridEntitySearch ( const GridPartType &gridPart )
      : hierarchicSearch_( gridPart.grid(), gridPart.indexSet() )
      {}

      EntityType operator() ( const GlobalCoordinateType &x ) const
      {
        return hierarchicSearch_.template findEntity< partition >( x );
      }

    private:
      Dune::HierarchicSearch< typename GridPartType::GridType, typename GridPartType::IndexSetType > hierarchicSearch_;
    };



    // EntitySearch
    // ------------

    template< class GridPart, int codim = 0, PartitionIteratorType partition = All_Partition >
    class EntitySearch
    : public std::conditional< GridPartCapabilities::hasGrid< GridPart >::v, GridEntitySearch< GridPart, codim, partition >, DefaultEntitySearch< GridPart, codim, partition > >::type
    {
      typedef EntitySearch< GridPart, codim, partition > ThisType;
      typedef typename std::conditional< GridPartCapabilities::hasGrid< GridPart >::v, GridEntitySearch< GridPart, codim, partition >, DefaultEntitySearch< GridPart, codim, partition > >::type BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit EntitySearch ( const GridPartType &gridPart )
      : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_ENTITYSEARCH_HH
