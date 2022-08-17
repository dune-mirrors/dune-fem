#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_HH

#include <cassert>

#include <dune/grid/common/gridview.hh>

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/extendedentity.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/geogridpart/capabilities.hh>
#include <dune/fem/gridpart/geogridpart/datahandle.hh>
#include <dune/fem/gridpart/geogridpart/entity.hh>
#include <dune/fem/gridpart/geogridpart/geometry.hh>
#include <dune/fem/gridpart/geogridpart/intersection.hh>
#include <dune/fem/gridpart/geogridpart/intersectioniterator.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>
#include <dune/fem/gridpart/idgridpart/iterator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class CoordFunction >
    class GeoGridPart;

    template< class CoordFunction >
    struct GeoGridPartFamily;



    // GeoGridPartData
    // ---------------

    template< class CoordFunction >
    struct GeoGridPartData
    {
      typedef CoordFunction CoordFunctionType;

      GeoGridPartData () = default;
      GeoGridPartData ( const CoordFunctionType &coordFunction ) : coordFunction_( &coordFunction ) {}

      operator const CoordFunctionType & () const { assert( coordFunction_ ); return *coordFunction_; }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
    };



    // GeoGridPartFamily
    // -----------------

    // Traits for dune-grid facades ("Gen-Gurke!")
    template< class CoordFunction >
    struct GeoGridPartFamily
    {
      typedef typename CoordFunction::RangeFieldType ctype;

      static const int dimension = CoordFunction::GridPartType::dimension;
      static const int dimensionworld = CoordFunction::FunctionSpaceType::dimRange;

      typedef GeoGridPartFamily< CoordFunction > GridPartFamily;

      struct Traits
      {
        typedef GeoGridPartData< CoordFunction > ExtraData;
        typedef CoordFunction CoordFunctionType;

        typedef typename CoordFunctionType::GridPartType HostGridPartType;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension - codim, dimensionworld, const GridPartFamily, GeoGeometry > Geometry;
          typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometry;

          typedef Dune::ExtendedEntity< codim, dimension, const GridPartFamily, GeoEntity > Entity;
          typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
        };

        typedef DeadIntersection< const GridPartFamily > IntersectionImplType;
        typedef DeadIntersectionIterator< const GridPartFamily > IntersectionIteratorImplType;

        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LeafIntersection;
        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LevelIntersection;

        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LeafIntersectionIterator;
        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const GridPartFamily, DeadIterator< typename Codim< 0 >::Entity > > HierarchicIterator;
      };

      template< int codim >
      struct Codim
      : public Traits::template Codim< codim >
      {};

      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      typedef typename Traits::HierarchicIterator HierarchicIterator;
    };



    template< class CoordFunction >
    struct GeoGridPartTraits
    {
      typedef GeoGridPart< CoordFunction > GridPartType;
      typedef GeoGridPartFamily< CoordFunction > GridPartFamily;
      typedef GeoGridPartFamily< CoordFunction > GridFamily;

      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::GridType GridType;

      //! type of twist utility
      typedef MetaTwistUtility< typename HostGridPartType :: TwistUtilityType >  TwistUtilityType;

      typedef IdIndexSet< const GridPartFamily > IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      typedef GeoIntersection< const GridPartFamily > IntersectionImplType;
      typedef GeoIntersectionIterator< const GridPartFamily > IntersectionIteratorImplType;

      typedef IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridPartFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridPartFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::Entity EntityType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridPartFamily, IdIterator< codim, pitype, const GridPartFamily > > IteratorType;
        };
      };

      typedef typename HostGridPartType::CommunicationType CommunicationType;

      static const bool conforming = HostGridPartType::Traits::conforming;
    };



    // GeoGridPart
    // -----------

    template< class CoordFunction >
    class GeoGridPart
    : public GridPartDefault< GeoGridPartTraits< CoordFunction > >
    {
      typedef GeoGridPart< CoordFunction > ThisType;
      typedef GridPartDefault< GeoGridPartTraits< CoordFunction > > BaseType;

      typedef typename GeoGridPartTraits< CoordFunction >::GridPartFamily GridPartFamily;

    public:
      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

      typedef CoordFunction CoordFunctionType;

      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::CommunicationType CommunicationType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit GeoGridPart ( const CoordFunctionType &coordFunction )
      : BaseType( const_cast< GridType& > (coordFunction.gridPart().grid() ) ),
        coordFunction_( &coordFunction ),
        indexSet_( hostGridPart().indexSet() )
      {}

      const IndexSetType &indexSet () const
      {
        return indexSet_;
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        return IdIterator< codim, pitype, const GridPartFamily >( coordFunction(), hostGridPart().template begin< codim, pitype >() );
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return IdIterator< codim, pitype, const GridPartFamily >( coordFunction(), hostGridPart().template end< codim, pitype >() );
      }

      int level () const
      {
        return hostGridPart().level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeoIntersectionIterator< const GridPartFamily >( entity, hostGridPart().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeoIntersectionIterator< const GridPartFamily >( entity, hostGridPart().iend( entity.impl().hostEntity() ) );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        GeoDataHandle< GridPartFamily, HostHandleType > handleWrapper( coordFunction(), handle );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      template< class LocalFunction >
      typename Codim< 0 >::EntityType
      exchangeGeometry ( const typename Codim< 0 >::EntityType &entity,
                         const LocalFunction &localCoordFunction ) const
      {
        return typename Codim< 0 >::EntityType::Implementation( entity.impl(), localCoordFunction );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return typename Codim< EntitySeed::codimension >::EntityType
                 ::Implementation( coordFunction(), hostGridPart().entity( seed ) );
      }

      // convert a grid entity to a grid part entity ("Gurke!")
      template< class Entity >
      MakeableInterfaceObject< typename Codim< Entity::codimension >::EntityType >
      convert ( const Entity &entity ) const
      {
        // make sure we have a grid entity
        const auto& gridEntity = Fem::gridEntity( entity );

        // create a grid part entity from a given grid entity
        typedef typename Codim< Entity::codimension >::EntityType EntityType;
        typedef typename EntityType::Implementation Implementation;
        typedef MakeableInterfaceObject< EntityType > EntityObj;

        // here, grid part information can be passed, if necessary
        return EntityObj( Implementation( coordFunction(), hostGridPart().convert( gridEntity ) ) );
      }

      // return reference to the coordfunction
      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_);
        return *coordFunction_;
      }

      // return reference to the host grid part
      const HostGridPartType &hostGridPart () const
      {
        return coordFunction().gridPart();
      }

      // return reference to the host grid part
      HostGridPartType &hostGridPart ()
      {
        return const_cast< HostGridPartType & >( coordFunction().gridPart() );
      }

    protected:
      const CoordFunctionType *coordFunction_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for GeoEntity
    // ------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::ExtendedEntity< codim, dim, GridFamily, GeoEntity > >
    {
      typedef Dune::ExtendedEntity< codim, dim, GridFamily, GeoEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };



    // EntitySearch for GeoGridPart
    // ----------------------------

    template< class CoordFunction, int codim, PartitionIteratorType partition >
    class EntitySearch< GeoGridPart< CoordFunction >, codim, partition >
    : public DefaultEntitySearch< GeoGridPart< CoordFunction >, codim, partition >
    {
      typedef EntitySearch< GeoGridPart< CoordFunction >, codim, partition > ThisType;
      typedef DefaultEntitySearch< GeoGridPart< CoordFunction >, codim, partition > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit EntitySearch ( const GridPartType &gridPart )
      : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH
