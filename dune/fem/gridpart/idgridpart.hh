#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_HH

#include <dune/grid/common/gridview.hh>

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/extendedentity.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/idgridpart/capabilities.hh>
#include <dune/fem/gridpart/idgridpart/datahandle.hh>
#include <dune/fem/gridpart/idgridpart/entity.hh>
#include <dune/fem/gridpart/idgridpart/geometry.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>
#include <dune/fem/gridpart/idgridpart/intersection.hh>
#include <dune/fem/gridpart/idgridpart/intersectioniterator.hh>
#include <dune/fem/gridpart/idgridpart/iterator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class HostGridPartImp >
    class IdGridPart;


    // IdGridPartTraits
    // ----------------

    template< class HostGridPartImp >
    struct IdGridPartTraits
    {
      typedef IdGridPart< HostGridPartImp > GridPartType;

      //! type of twist utility
      typedef MetaTwistUtility< typename HostGridPartImp :: TwistUtilityType >  TwistUtilityType;

      // Traits for dune-grid facades ("Gen-Gurke!")
      struct GridFamily
      {
        typedef typename HostGridPartImp::ctype ctype;

        static const int dimension = HostGridPartImp::dimension;
        static const int dimensionworld = HostGridPartImp::dimensionworld;

        struct Traits
        {
          typedef HostGridPartImp HostGridPartType;

          struct EmptyData {};

          // type of data passed to entities, intersections, and iterators
          // for IdGridPart this is just an empty place holder
          typedef EmptyData ExtraData;

          template< int codim >
          struct Codim
          {
            typedef Dune::Geometry< dimension - codim, dimensionworld, const GridFamily, IdGeometry > Geometry;
            typedef Dune::Geometry< dimension - codim, dimension, const GridFamily, IdLocalGeometry > LocalGeometry;

            typedef Dune::ExtendedEntity< codim, dimension, const GridFamily, IdEntity > Entity;
            typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
          };

          typedef DeadIntersection< const GridFamily > IntersectionImplType;
          typedef DeadIntersectionIterator< const GridFamily > IntersectionIteratorImplType;

          typedef Dune::Intersection< const GridFamily, IntersectionImplType > LeafIntersection;
          typedef Dune::Intersection< const GridFamily, IntersectionImplType > LevelIntersection;

          typedef Dune::IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > LeafIntersectionIterator;
          typedef Dune::IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > LevelIntersectionIterator;

          typedef Dune::EntityIterator< 0, const GridFamily, DeadIterator< typename Codim< 0 >::Entity > > HierarchicIterator;
        };

        template< int codim >
        struct Codim
        : public Traits::template Codim< codim >
        {};

        typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
        typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

        typedef typename Traits::HierarchicIterator HierarchicIterator;
      };
      typedef typename GridFamily::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::GridType GridType;

      typedef IdIndexSet< const GridFamily > IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      typedef IdIntersectionIterator < const GridFamily > IntersectionIteratorImplType;
      typedef IdIntersection< const GridFamily > IntersectionImplType;
      typedef IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridFamily::Traits::template Codim< codim >::Entity EntityType;
        typedef typename GridFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridFamily, IdIterator< codim, pitype, const GridFamily > > IteratorType;
        };
      };

      typedef typename HostGridPartType::CommunicationType CommunicationType;

      static const bool conforming = HostGridPartType::Traits::conforming;
    };



    // IdGridPart
    // ----------

    template< class HostGridPartImp >
    class IdGridPart
    : public GridPartDefault< IdGridPartTraits< HostGridPartImp > >
    {
      typedef IdGridPart< HostGridPartImp > ThisType;
      typedef GridPartDefault< IdGridPartTraits< HostGridPartImp > > BaseType;

      typedef typename IdGridPartTraits< HostGridPartImp >::GridFamily GridFamily;

    public:
      typedef typename GridFamily::Traits::HostGridPartType HostGridPartType;

      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::CommunicationType CommunicationType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit IdGridPart ( GridType &grid )
      : BaseType( grid ),
        hostGridPart_( grid ),
        indexSet_( hostGridPart_.indexSet() )
      {}

      explicit IdGridPart ( const IdGridPart &other )
      : BaseType( other ),
        hostGridPart_( other.hostGridPart() ),
        indexSet_( hostGridPart_.indexSet() )
      {}

      explicit IdGridPart ( const HostGridPartType &hostGridPart )
      : BaseType( const_cast< GridType& > ( hostGridPart.grid() ) ),
        hostGridPart_( hostGridPart ),
        indexSet_( hostGridPart_.indexSet() )
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
        return IdIterator< codim, pitype, const GridFamily >( data(), hostGridPart().template begin< codim, pitype >() );
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
        return IdIterator< codim, pitype, const GridFamily >( data(), hostGridPart().template end< codim, pitype >() );
      }

      int level () const
      {
        return hostGridPart().level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return IdIntersectionIterator< const GridFamily >( data(), hostGridPart().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return IdIntersectionIterator< const GridFamily >( data(), hostGridPart().iend( entity.impl().hostEntity() ) );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        IdDataHandle< HostHandleType, GridFamily > handleWrapper( data(), handle );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        typedef typename Codim< EntitySeed::codimension >::EntityType::Implementation EntityImp;
        return EntityImp( data(), hostGridPart().entity( seed ) );
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
        return EntityObj( Implementation( data(), hostGridPart().convert( gridEntity ) ) );
      }

      const HostGridPartType &hostGridPart () const { return hostGridPart_; }

      HostGridPartType &hostGridPart () { return hostGridPart_; }

      typedef typename GridFamily::Traits::ExtraData ExtraData;
      ExtraData data () const { return ExtraData(); }

    protected:
      HostGridPartType hostGridPart_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for IdEntity
    // -----------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::ExtendedEntity< codim, dim, GridFamily, IdEntity > >
    {
      typedef Dune::ExtendedEntity< codim, dim, GridFamily, IdEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };



    // EntitySearch for IdGridPart
    // ---------------------------

    template< class HostGridPart, int codim, PartitionIteratorType partition >
    class EntitySearch< IdGridPart< HostGridPart >, codim, partition >
    {
      typedef EntitySearch< IdGridPart< HostGridPart >, codim, partition > ThisType;

    public:
      typedef IdGridPart< HostGridPart > GridPartType;
      typedef typename GridPartType::ExtraData  ExtraData;

      typedef typename GridPartType::template Codim< codim >::EntityType EntityType;

      typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinateType;

      explicit EntitySearch ( const GridPartType &gridPart )
      : hostEntitySearch_( gridPart.hostGridPart() ),
        data_( gridPart.data() )
      {}

      EntityType operator() ( const GlobalCoordinateType &x ) const
      {
        typedef typename EntityType::Implementation EntityImpl;
        return EntityImpl( data_, hostEntitySearch_( x ) );
      }

    protected:
      const EntitySearch< HostGridPart > hostEntitySearch_;
      ExtraData data_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_HH
