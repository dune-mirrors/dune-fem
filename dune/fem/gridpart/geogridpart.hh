#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_HH

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/geogridpart/capabilities.hh>
#include <dune/fem/gridpart/geogridpart/datahandle.hh>
#include <dune/fem/gridpart/geogridpart/entity.hh>
#include <dune/fem/gridpart/geogridpart/entitypointer.hh>
#include <dune/fem/gridpart/geogridpart/intersection.hh>
#include <dune/fem/gridpart/geogridpart/intersectioniterator.hh>
#include <dune/fem/gridpart/geogridpart/iterator.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class CoordFunction >
    class GeoGridPart;

    template< class CoordFunction >
    class GeoGridPartFamily;

  } // namespace Fem



  namespace GenericGeometry
  {

    // GeometryTraits for GeoGridPart
    // ------------------------------
  
    template< class CoordFunction >
    struct GlobalGeometryTraits< Fem::GeoGridPartFamily< CoordFunction > >
    : public DefaultGeometryTraits< typename CoordFunction::RangeFieldType, CoordFunction::FunctionSpaceType::dimDomain, CoordFunction::FunctionSpaceType::dimRange >
    {
      typedef Fem::GeoGridPartFamily< CoordFunction > GridPartFamily;

      typedef DuneCoordTraits< typename CoordFunction::RangeFieldType > CoordTraits;

      static const int dimGrid = CoordFunction::FunctionSpaceType::dimDomain;
      static const int dimWorld = CoordFunction::FunctionSpaceType::dimRange;

      static const bool hybrid = !Capabilities::hasSingleGeometryType< typename CoordFunction::GridType >::v;
      // this value is only used when hybrid is false (and only valid in that case)
      static const unsigned int topologyId = Capabilities::hasSingleGeometryType< typename CoordFunction::GridType >::topologyId;

      template< class Topology >
      struct Mapping
      {
        typedef Fem::GeoCornerStorage< Topology, const GridPartFamily > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage > type;
      };

      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
        static const EvaluationType evaluateNormal = ComputeOnDemand;
      };
    };
    
  } // namespace GenericGeometry



  namespace FacadeOptions
  {

    template< int mydim, int cdim, class CoordFunction >
    struct StoreGeometryReference< mydim, cdim, const Fem::GeoGridPartFamily< CoordFunction >, Dune::GenericGeometry::Geometry >
    {
      static const bool v = false;
    };

  } // namespace FacadeOptions



  namespace Fem
  {

    // GeoGridPartFamily
    // -----------------

    // Traits for dune-grid facades ("Gen-Gurke!")
    template< class CoordFunction >
    struct GeoGridPartFamily
    {
      typedef typename CoordFunction::RangeFieldType ctype;

      static const int dimension = CoordFunction::FunctionSpaceType::dimDomain;
      static const int dimensionworld = CoordFunction::FunctionSpaceType::dimRange;

      typedef GeoGridPartFamily< CoordFunction > GridPartFamily;

      struct Traits
      {
        typedef CoordFunction CoordFunctionType;
        typedef typename CoordFunctionType::GridPartType HostGridPartType;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension - codim, dimensionworld, const GridPartFamily, Dune::GenericGeometry::Geometry > Geometry;
          typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometry;

          typedef GeoEntityPointer< GeoEntityPointerTraits< codim, const GridPartFamily > > EntityPointerImpl;
          typedef Dune::EntityPointer< const GridPartFamily, EntityPointerImpl > EntityPointer;

          typedef Dune::Entity< codim, dimension, const GridPartFamily, GeoEntity > Entity;
          typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
        };

        typedef Dune::Intersection< const GridPartFamily, DeadIntersection > LeafIntersection;
        typedef Dune::Intersection< const GridPartFamily, DeadIntersection > LevelIntersection;

        typedef Dune::IntersectionIterator< const GridPartFamily, DeadIntersectionIterator, DeadIntersection > LeafIntersectionIterator;
        typedef Dune::IntersectionIterator< const GridPartFamily, DeadIntersectionIterator, DeadIntersection > LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const GridPartFamily, DeadIterator< typename Codim< 0 >::EntityPointerImpl > > HierarchicIterator;
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

      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::GridType GridType;

      typedef IdIndexSet< typename HostGridPartType::IndexSetType > IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      typedef IntersectionIterator< const GridPartFamily, GeoIntersectionIterator, GeoIntersection > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridPartFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridPartFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::EntityPointer EntityPointerType;
        typedef typename GridPartFamily::Traits::template Codim< codim >::Entity EntityType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridPartFamily, GeoIterator< codim, pitype, const GridPartFamily > > IteratorType;
        };
      };

      static const bool conforming = HostGridPartType::conforming;
    };



    // GeoGridPart
    // -----------

    template< class CoordFunction >
    class GeoGridPart
    : public GridPartInterface< GeoGridPartTraits< CoordFunction > >
    {
      typedef GeoGridPart< CoordFunction > ThisType;
      typedef GridPartInterface< GeoGridPartTraits< CoordFunction > > BaseType;

      typedef typename GeoGridPartTraits< CoordFunction >::GridPartFamily GridPartFamily;

      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

    public:
      typedef CoordFunction CoordFunctionType;

      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit GeoGridPart ( const CoordFunctionType &coordFunction )
      : coordFunction_( coordFunction ),
        indexSet_( hostGridPart().indexSet() )
      {}

      const GridType &grid () const
      {
        return hostGridPart().grid();
      }

      GridType &grid ()
      { 
        return const_cast< GridType & >( hostGridPart().grid() );
      }

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
        return GeoIterator< codim, pitype, const GridPartFamily >( coordFunction_, hostGridPart().template begin< codim, pitype >() );
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
        return GeoIterator< codim, pitype, const GridPartFamily >( coordFunction_, hostGridPart().template end< codim, pitype >() );
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

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return hostGridPart().boundaryId( intersection.impl().hostIntersection() );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        GeoDataHandle< GridPartFamily, HostHandleType > handleWrapper( coordFunction_, handle );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      template< class LCFTraits >
      typename Codim< 0 >::EntityPointerType
      exchangeGeometry ( const typename Codim< 0 >::EntityType &entity,
                         const LocalFunction< LCFTraits > &localCoordFunction ) const
      {
        return typename Codim< 0 >::EntityPointerType::Implementation( entity.impl(), localCoordFunction );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityPointerType
      entityPointer ( const EntitySeed &seed ) const
      {
        return Codim< EntitySeed::codimension >::EntityPointerType
                 ::Implementation( coordFunction_, hostGridPart().entityPointer( seed ) );
      }

      // convert a grid entity to a grid part entity ("Gurke!")
      template< class Entity > 
      MakeableInterfaceObject< typename Codim< Entity::codimension >::EntityType >
      convert ( const Entity &entity ) const
      {
        // create a grid part entity from a given grid entity
        typedef typename Codim< Entity::codimension >::EntityType EntityType;
        typedef typename EntityType::Implementation Implementation;
        typedef MakeableInterfaceObject< EntityType > EntityObj;

        // here, grid part information can be passed, if necessary
        return EntityObj( Implementation( coordFunction_, entity ) );
      }

    private:
      const HostGridPartType &hostGridPart () const
      {
        return coordFunction_.gridPart();
      }

      const CoordFunctionType &coordFunction_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for GeoEntity
    // ------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::Entity< codim, dim, GridFamily, GeoEntity > >
    {
      typedef Dune::Entity< codim, dim, GridFamily, GeoEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH
