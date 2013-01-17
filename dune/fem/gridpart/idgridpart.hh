#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_HH

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/idgridpart/capabilities.hh>
#include <dune/fem/gridpart/idgridpart/datahandle.hh>
#include <dune/fem/gridpart/idgridpart/entity.hh>
#include <dune/fem/gridpart/idgridpart/entitypointer.hh>
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

    template< class HostGridPart >
    class IdGridPart;


    // IdGridPartTraits
    // ----------------

    template< class HostGridPart >
    struct IdGridPartTraits
    {
      typedef IdGridPart< HostGridPart > GridPartType;

      //! type of twist utility 
      typedef MetaTwistUtility< typename HostGridPart :: TwistUtilityType >  TwistUtilityType; 

      // Traits for dune-grid facades ("Gen-Gurke!")
      struct GridFamily
      {
        typedef typename HostGridPart::ctype ctype;

        static const int dimension = HostGridPart::dimension;
        static const int dimensionworld = HostGridPart::dimensionworld;

        struct Traits
        {
          typedef HostGridPart HostGridPartType;

          template< int codim >
          struct Codim
          {
            typedef Dune::Geometry< dimension - codim, dimensionworld, const GridFamily, IdGeometry > Geometry;
            typedef Dune::Geometry< dimension - codim, dimension, const GridFamily, IdLocalGeometry > LocalGeometry;

            typedef IdEntityPointer< IdEntityPointerTraits< codim, const GridFamily > > EntityPointerImpl;
            typedef Dune::EntityPointer< const GridFamily, EntityPointerImpl > EntityPointer;

            typedef Dune::Entity< codim, dimension, const GridFamily, IdEntity > Entity;
            typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
          };

          typedef DeadIntersection< const GridFamily > IntersectionImplType;
          typedef DeadIntersectionIterator< const GridFamily > IntersectionIteratorImplType;

          typedef Dune::Intersection< const GridFamily, IntersectionImplType > LeafIntersection;
          typedef Dune::Intersection< const GridFamily, IntersectionImplType > LevelIntersection;

          typedef Dune::IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > LeafIntersectionIterator;
          typedef Dune::IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > LevelIntersectionIterator;

          typedef Dune::EntityIterator< 0, const GridFamily, DeadIterator< typename Codim< 0 >::EntityPointerImpl > > HierarchicIterator;
        };

        template< int codim >
        struct Codim
        : public Traits::template Codim< codim >
        {};

        typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
        typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

        typedef typename Traits::HierarchicIterator HierarchicIterator;
      };

      typedef typename HostGridPart::GridType GridType;

      typedef IdIndexSet< typename HostGridPart::IndexSetType > IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPart::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPart::indexSetInterfaceType;

      typedef IdIntersectionIterator < const GridFamily > IntersectionIteratorImplType;
      typedef IdIntersection< const GridFamily > IntersectionImplType;
      typedef IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridFamily::Traits::template Codim< codim >::EntityPointer EntityPointerType;
        typedef typename GridFamily::Traits::template Codim< codim >::Entity EntityType;

        typedef typename GridFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridFamily, IdIterator< codim, pitype, const GridFamily > > IteratorType;
        };
      };

      static const bool conforming = HostGridPart::conforming;
    };



    // IdGridPart
    // ----------

    template< class HostGridPart >
    class IdGridPart
    : public GridPartInterface< IdGridPartTraits< HostGridPart > >
    {
      typedef IdGridPart< HostGridPart > ThisType;
      typedef GridPartInterface< IdGridPartTraits< HostGridPart > > BaseType;

      typedef typename IdGridPartTraits< HostGridPart >::GridFamily GridFamily;

    public:
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit IdGridPart ( GridType &grid )
      : hostGridPart_( grid ),
        indexSet_( hostGridPart_.indexSet() )
      {}

      explicit IdGridPart ( const HostGridPart &hostGridPart )
      : hostGridPart_( hostGridPart ),
        indexSet_( hostGridPart_.indexSet() )
      {}

      const GridType &grid () const
      {
        return hostGridPart_.grid();
      }

      GridType &grid ()
      { 
        return hostGridPart_.grid();
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
        return IdIterator< codim, pitype, const GridFamily >( hostGridPart_.template begin< codim, pitype >() );
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
        return IdIterator< codim, pitype, const GridFamily >( hostGridPart_.template end< codim, pitype >() );
      }

      int level () const 
      {
        return hostGridPart_.level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return IdIntersectionIterator< const GridFamily >( hostGridPart_.ibegin( entity.impl().hostEntity() ) );
      }
      
      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const 
      {
        return IdIntersectionIterator< const GridFamily >( hostGridPart_.iend( entity.impl().hostEntity() ) );
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return hostGridPart_.boundaryId( intersection.impl().hostIntersection() );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        IdDataHandle< GridFamily, HostHandleType > handleWrapper( handle );
        hostGridPart_.communicate( handleWrapper, iftype, dir );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityPointerType
      entityPointer ( const EntitySeed &seed ) const
      {
        typedef typename Codim< EntitySeed::codimension >::EntityPointerType::Implementation EntityPointerImp;
        return EntityPointerImp( hostGridPart_.entityPointer( seed ) );
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
        return EntityObj( Implementation( entity ) );
      }

    private:
      HostGridPart hostGridPart_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for IdEntity
    // -----------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::Entity< codim, dim, GridFamily, IdEntity > >
    {
      typedef Dune::Entity< codim, dim, GridFamily, IdEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_HH
