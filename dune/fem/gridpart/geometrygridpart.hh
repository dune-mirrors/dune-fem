#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH

#include <cassert>

#include <dune/common/version.hh>

#include <dune/fem/function/localfunction/const.hh>

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/extendedentity.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>
#include <dune/fem/gridpart/idgridpart/iterator.hh>

#include <dune/fem/gridpart/common/compositegeometry.hh>
#include <dune/fem/gridpart/common/localfunctiongeometry.hh>
#include <dune/fem/gridpart/common/sharedgeometry.hh>
#include <dune/fem/gridpart/geometrygridpart/capabilities.hh>
#include <dune/fem/gridpart/geometrygridpart/entity.hh>
#include <dune/fem/gridpart/geometrygridpart/datahandle.hh>
#include <dune/fem/gridpart/geometrygridpart/intersection.hh>
#include <dune/fem/gridpart/geometrygridpart/intersectioniterator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridFunctionType >
    class GeometryGridPart;



    // GeometryGridPartData
    // --------------------

    template< class GridFunction >
    struct GeometryGridPartData
    {
      typedef GridFunction GridFunctionType;

      GeometryGridPartData () noexcept = default;
      GeometryGridPartData ( const GridFunctionType &gridFunction ) noexcept : gridFunction_( &gridFunction ) {}

      operator const GridFunctionType & () const { assert( gridFunction_ ); return *gridFunction_; }

    private:
      const GridFunctionType *gridFunction_ = nullptr;
    };



    // GeometryGridPartFamily
    // ----------------------

    template< class GridFunction >
    struct GeometryGridPartFamily
    {
      typedef GridFunction GridFunctionType;
      typedef typename GridFunction::RangeFieldType ctype;

      static const int dimension = GridFunction::GridPartType::dimension;
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      typedef GeometryGridPartFamily< GridFunction > GridPartFamily;

      struct Traits
      {
        typedef GeometryGridPartData< GridFunction > ExtraData;
        typedef GridFunction GridFunctionType;
        typedef typename GridFunctionType::GridPartType HostGridPartType;

        typedef SharedGeometry< LocalFunctionGeometry< ConstLocalFunction<GridFunction> > > ElementGeometryImpl;

        template< int codim >
        struct Codim
        {
          typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometry;

          template< int mydim, int cdim, class Grid >
          using GeometryImpl = std::conditional_t< mydim == dimension, ElementGeometryImpl, CompositeGeometry< ElementGeometryImpl, LocalGeometry > >;

          typedef Dune::Geometry< dimension - codim, dimensionworld, const GridPartFamily, GeometryImpl > Geometry;

          typedef Dune::ExtendedEntity< codim, dimension, const GridPartFamily, GeometryGridPartEntity > Entity;
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



    // GeometryGridPartTraits
    // ----------------------

    template< class GridFunction >
    struct GeometryGridPartTraits
    {
      typedef GridFunction GridFunctionType;
      typedef typename GridFunction::GridPartType HostGridPartType;
      typedef GeometryGridPart< GridFunction > GridPartType;
      typedef GeometryGridPartFamily< GridFunction > GridPartFamily;
      typedef GeometryGridPartFamily< GridFunction > GridFamily;

      typedef GridPartType GridViewType;

      static const int dimension = GridFunction::GridPartType::dimension;
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      //! type of twist utility
      typedef MetaTwistUtility< typename HostGridPartType::TwistUtilityType >  TwistUtilityType;

      typedef IdIndexSet< const GridPartFamily > IndexSetType;

      typedef typename HostGridPartType::GridType GridType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      typedef GeometryGridPartIntersectionIterator< const GridFamily > IntersectionIteratorImplType;
      typedef GeometryGridPartIntersection< const GridFamily > IntersectionImplType;
      typedef IntersectionIterator< const GridFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim : public GridFamily::Traits::template Codim< codim >
      {
        typedef typename GridFamily::Traits::template Codim< codim > BaseType;

        typedef typename BaseType::Geometry       GeometryType;
        typedef typename BaseType::LocalGeometry  LocalGeometryType;

        typedef typename BaseType::Entity         EntityType;
        typedef typename BaseType::EntitySeed     EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridFamily, IdIterator< codim, pitype, const GridFamily > > IteratorType;
          typedef IteratorType Iterator;
        };
      };

      typedef typename HostGridPartType::CommunicationType CommunicationType;
      typedef CommunicationType Communication;
      static const bool conforming = HostGridPartType::Traits::conforming;
    };



    // GeometryGridPart
    // ----------------

    template< class GridFunction >
    class GeometryGridPart
      : public GridPartDefault< GeometryGridPartTraits< GridFunction > >
    {
    public:
      typedef GridFunction GridFunctionType;

    private:
      typedef GeometryGridPart< GridFunctionType > ThisType;
      typedef GridPartDefault< GeometryGridPartTraits< GridFunctionType > > BaseType;
      typedef typename GeometryGridPartTraits< GridFunctionType >::GridFamily GridFamily;

    public:
      typedef typename GridFunctionType::GridPartType HostGridPartType;
      //! \brief type of grid
      typedef typename BaseType::GridType GridType;
      //! \brief type of grid
      typedef typename BaseType::Grid     Grid;
      //! \brief index set use in this gridpart
      typedef typename BaseType::IndexSetType IndexSetType;
      //! \brief index set use in this gridpart
      typedef typename BaseType::IndexSet     IndexSet;
      //! \brief type of intersection iterator
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      //! \brief type of intersection iterator
      typedef typename BaseType::IntersectionIterator     IntersectionIterator;
      //! \brief type of intersection
      typedef typename BaseType::IntersectionType IntersectionType;
      //! \brief type of intersection
      typedef typename BaseType::Intersection Intersection;
      //! \brief Collective communication
      typedef typename BaseType::CommunicationType CommunicationType;
      //! \brief Collective communication
      typedef typename BaseType::Communication     Communication;
      typedef typename BaseType::GridViewType GridViewType;

      // the interface takes this from the grid
      static const int dimensionworld = GridFunction::FunctionSpaceType::dimRange;

      template< int codim >
      struct Codim
        : public BaseType::template Codim< codim >
      {};

      explicit GeometryGridPart ( const GridFunctionType &gridFunction )
        : BaseType( const_cast< GridType& > (gridFunction.gridPart().grid()) ),
          gridFunction_( &gridFunction ),
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
        return IdIterator< codim, pitype, const GridFamily >( gridFunction(), hostGridPart().template begin< codim, pitype >() );
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
        return IdIterator< codim, pitype, const GridFamily >( gridFunction(), hostGridPart().template end< codim, pitype >() );
      }

      int level () const
      {
        return hostGridPart().level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeometryGridPartIntersectionIterator< const GridFamily >( entity, hostGridPart().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeometryGridPartIntersectionIterator< const GridFamily >( entity, hostGridPart().iend( entity.impl().hostEntity() ) );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        GeometryGridPartDataHandle< GridFamily, HostHandleType > handleWrapper( handle, gridFunction() );
        hostGridPart().communicate( handleWrapper, iftype, dir );
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
        return EntityObj( Implementation( gridFunction(), hostGridPart().convert( gridEntity ) ) );
      }
      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return convert( hostGridPart().entity(seed) );
      }

      const HostGridPartType &hostGridPart () const
      {
        return gridFunction().gridPart();
      }

      const GridFunctionType& gridFunction() const
      {
        assert( gridFunction_ );
        return *gridFunction_;
      }

    protected:
      const GridFunctionType *gridFunction_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for GeometryGridPartEntity
    // -------------------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::ExtendedEntity< codim, dim, GridFamily, GeometryGridPartEntity > >
    {
      typedef Dune::ExtendedEntity< codim, dim, GridFamily, GeometryGridPartEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };



    // EntitySearch for GeometryGridPart
    // ---------------------------------

    template< class GridFunction, int codim, PartitionIteratorType partition >
    class EntitySearch< GeometryGridPart< GridFunction >, codim, partition >
      : public DefaultEntitySearch< GeometryGridPart< GridFunction >, codim, partition >
    {
      typedef EntitySearch< GeometryGridPart< GridFunction >, codim, partition > ThisType;
      typedef DefaultEntitySearch< GeometryGridPart< GridFunction >, codim, partition > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit EntitySearch ( const GridPartType &gridPart )
        : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_HH
