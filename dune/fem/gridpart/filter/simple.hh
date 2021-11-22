#ifndef DUNE_FEM_GRIDPART_FILTER_SIMPLE_HH
#define DUNE_FEM_GRIDPART_FILTER_SIMPLE_HH

#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/fem/gridpart/filter/filter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridPart >
    struct SimpleFilter;



    // SimpleFilterTraits
    // ---------------------

    template< class GridPart >
    struct SimpleFilterTraits
    {
      typedef SimpleFilter< GridPart > FilterType;

      template< int codim >
      struct Codim
      {
        typedef typename GridPart::template Codim< codim >::EntityType EntityType;
      };
    };



    // SimpleFilter
    // ---------------

    template< class GridPart >
    class SimpleFilter
      : public FilterDefaultImplementation< SimpleFilterTraits< GridPart > >
    {
      typedef SimpleFilter< GridPart > This;
      typedef FilterDefaultImplementation< SimpleFilterTraits< GridPart > > Base;

      typedef std::make_integer_sequence< int, GridPart::dimension+1 > AllCodims;

    public:
      template< int codim >
      using Codim = typename Base::template Codim< codim >;

      template< class Contains >
      SimpleFilter ( const GridPart &gridPart, Contains contains, int domainId )
        : mapper_( gridPart, [](Dune::GeometryType,int ) {return true;} ),
          contains_( mapper_.size(), false ),
          mapper0_( gridPart, [](Dune::GeometryType gt,int dim) {return gt.dim()==dim;} ),
          domainIds_( mapper0_.size(), -1 )
      {
        for( const typename Codim< 0 >::EntityType &entity : elements( gridPart ) )
        {
          domainIds_[ mapper0_.index( entity ) ] = contains( entity );
          if( domainIds_[ mapper0_.index( entity ) ] == domainId )
            mark( entity, AllCodims() );
        }
      }

      template< int codim >
      bool contains ( const typename Codim< codim >::EntityType &entity ) const
      {
        return contains_[ mapper_.index( entity ) ];
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return contains< Entity::codimension >( entity );
      }

      template< class Intersection >
      bool interiorIntersection ( const Intersection &intersection ) const
      {
        return contains( intersection.outside() );
      }

      template< class Intersection >
      bool intersectionBoundary ( const Intersection & ) const
      {
        return true;
      }

      template< class Intersection >
      int intersectionBoundaryId ( const Intersection &intersection ) const
      {
        return ( intersection.boundary() ?
            Dune::Fem::BoundaryIdProvider<typename GridPart::GridType>::boundaryId( intersection ) :
                 domainIds_[ mapper0_.index( intersection.outside() ) ] );
      }

      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & ) const
      {
        return false;
      }

    private:
      template< int codim >
      void mark ( const typename Codim< 0 >::EntityType &entity, Dune::Codim< codim > )
      {
        for( unsigned int i = 0, n = entity.subEntities( codim ); i < n; ++i )
          contains_[ mapper_.subIndex( entity, i, codim ) ] = true;
      }

      template< int... codim >
      void mark ( const typename Codim< 0 >::EntityType &entity, std::integer_sequence< int, codim... > )
      {
        std::ignore = std::make_tuple( (mark( entity, Dune::Codim< codim >() ), codim)... );
      }

      MultipleCodimMultipleGeomTypeMapper< typename GridPart::GridViewType > mapper_;
      std::vector< bool > contains_;
      MultipleCodimMultipleGeomTypeMapper< typename GridPart::GridViewType > mapper0_;
      std::vector< int > domainIds_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_SIMPLE_HH
