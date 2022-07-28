#ifndef DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH
#define DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH

//- dune-common includes
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/exceptions.hh>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/grid.hh>

//- dune-fem includes
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/misc/boundaryidprovider.hh>

namespace Dune
{

  namespace Fem
  {

    /** \addtogroup GridPart
     *
     * Grid parts allow to define a view on a given DUNE grid, treating the
     * underlying grid as a container for entities.
     *
     * All parts of the dune-fem package rely on grid parts to access the entities
     * of the grid. For example, discrete functions are defined on the set of
     * entities accesseable by the given GridPart implementation using the iterator
     * and index set provided by the GridPart.
     *
     * \section GridPart Interface and available Implementations
     *
     * The interface for a GridPart is implemented by the class template
     * GridPartInterface. Basically, a GridPart provides the following
     * functionality:
     * - The underlying grid can be accessed through the grid method.
     * - The indexSet method provides a suitable dune-fem index set for the grid
     *   part.
     * - Pairs of begin / end methods provide iterators over the entities of a
     *   given codimension belonging to the grid part.
     * - A pair of ibegin / iend methods provide suitable intersection iterators
     *   for a given entity of codimension 0.
     * - For parallel computations, a suitable communicate method is provided.
     * .
     *
     * The following grid parts have been implemented:
     * - LeafGridPart: A view of the leaf grid,
     * - LevelGridPart: A view of a given grid level,
     * - FilteredGridPart: A view filtering another grid part.
     *
     * \todo Implement a grid part for a given grid view (Suggestion: use the
     *       name GridPart).
     */



    /**
     * @addtogroup GridPart
     *
     * @{
     */

    //! \brief Interface for the GridPart classes
    //! A GridPart class allows to access only a specific subset of a grid's
    //! entities. A GridPart implementation provides the corresponding index set
    //! and a begin/end iterator pair for accessing those entities, the
    //! corresponding intersection iterators and a appropriate communication
    //! method.
    //! GridParts are used to parametrize spaces (see DiscreteFunctionSpaceDefault [in dune-fem]).
    template< class GridPartTraits >
    class GridPartInterface
    {
      typedef GridPartInterface< GridPartTraits > ThisType;

    public:
      //! \brief Type of the Traits
      typedef GridPartTraits Traits;

      //! \brief Type of the implementation
      typedef typename Traits::GridPartType GridPartType;

      //! \brief type of Grid implementation
      typedef typename Traits::GridType GridType;
      //! \brief type of Grid implementation
      typedef GridType Grid;

      //! \brief Index set implementation
      typedef typename Traits::IndexSetType IndexSetType;
      //! \brief Index set implementation
      typedef IndexSetType IndexSet;

      //! \brief Collective communication
      typedef typename Traits::CommunicationType CommunicationType;
      //! \brief Collective communication
      typedef CommunicationType Communication;

      //! \brief Twist utility type
      typedef typename Traits::TwistUtilityType TwistUtilityType;

      //! \brief Maximum Partition type, the index set provides indices for
      static const PartitionIteratorType indexSetPartitionType
        = Traits::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType
        = Traits::indexSetInterfaceType;

      //! \brief type of IntersectionIterator
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
      //! \brief type of IntersectionIterator
      typedef IntersectionIteratorType IntersectionIterator;

      //! \brief type of Intersection
      typedef typename IntersectionIteratorType::Intersection IntersectionType;
      //! \brief type of Intersection
      typedef IntersectionType Intersection;

      typedef GridPartType GridViewType;

      typedef typename GridType::ctype ctype;

      static const int dimension = GridType::dimension;
      static const int dimensionworld = GridType::dimensionworld;

      static const bool conforming = Traits::conforming ;

      template< int codim >
      struct Codim
      {
        typedef typename Traits::template Codim< codim >::GeometryType       GeometryType;
        typedef typename Traits::template Codim< codim >::LocalGeometryType  LocalGeometryType;

        typedef typename Traits::template Codim< codim >::EntityType         EntityType;
        typedef typename Traits::template Codim< codim >::EntitySeedType     EntitySeedType;

        // GridView typedefs interface
        typedef GeometryType       Geometry;
        typedef LocalGeometryType  LocalGeometry;
        typedef EntityType         Entity;
        typedef EntitySeedType     EntitySeed;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
            IteratorType;
          // GridView typedef
          typedef IteratorType Iterator;
        };

        typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;
        // GridView typedef
        typedef IteratorType Iterator;
      };

      //! \brief Returns const reference to the underlying grid
      const GridType &grid () const
      {
        CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
        return asImp().grid();
      }
      //! \brief Returns reference to the underlying grid
      GridType &grid ()
      {
        CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
        return asImp().grid();
      }

      //! \brief Returns reference to index set of the underlying grid
      const IndexSetType& indexSet() const
      {
        CHECK_INTERFACE_IMPLEMENTATION((asImp().indexSet()));
        return asImp().indexSet();
      }

      //! \brief returns true if current grid part represents a conforming grid
      bool isConforming() const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().isConforming()) );
        return asImp().isConforming();
      }

      /** \brief obtain number of entities in a given codimension */
      int size ( int codim ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().size( codim )) );
        return asImp().size( codim );
      }

      /** \brief obtain number of entities with a given geometry type */
      int size ( const GeometryType &type ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().size( type )) );
        return asImp().size( type );
      }

      /** \brief obtain begin iterator for the interior-border partition
       *
       *  \tparam  codim  codimension for which the iterator is requested
       */
      template< int codim >
      typename Codim< codim >::IteratorType
      begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim >()) );
        return asImp().template begin< codim >();
      }

      /** \brief obtain begin iterator for the given partition
       *
       *  \tparam  codim   codimension for which the iterator is requested
       *  \tparam  pitype  requested partition iterator type
       */
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim, pitype >()) );
        return asImp().template begin< codim, pitype >();
      }

      /** \brief obtain end iterator for the interior-border partition
       *
       *  \tparam  codim  codimension for which the iterator is requested
       */
      template< int codim >
      typename Codim< codim >::IteratorType
      end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim >()) );
        return asImp().template end< codim >();
      }

      /** \brief obtain end iterator for the given partition
       *
       *  \tparam  codim   codimension for which the iterator is requested
       *  \tparam  pitype  requested partition iterator type
       */
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim, pitype >()) );
        return asImp().template end< codim, pitype >();
      }

      //! \brief Level of the grid part
      [[deprecated("Do not use this since it is no longer part of the interface")]]
      int level () const
      {
        DUNE_THROW(NotImplemented,"GridPart::level has been removed!");
        return -1;
      }

      //! \brief ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType
      ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().ibegin( entity )) );
        return asImp().ibegin( entity );
      }

      //! \brief iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( (asImp().iend( entity )) );
        return asImp().iend( entity );
      }

      //! \brief return boundary if given an intersection
      [[deprecated("Use BoundaryIdProvider::boundaryId directly!")]]
      int boundaryId ( const IntersectionType &intersection ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().boundaryId( intersection ) );
        return asImp().boundaryId( intersection );
      }

      /** \brief obtain collective communication object */
      const CommunicationType &comm () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().comm() );
        return asImp().comm();
      }

      //! \brief corresponding communication method for grid part
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( (asImp().communicate( data, iftype, dir )) );
      }

      //! \brief obtain entity pointer from entity seed
      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().entity( seed ) );
        return asImp().entity( seed );
      }

      /*! \brief convert the grid's entity to a grid part entity
          Usually the parameter is GridType :: Codim< codim > :: Entity
          and the return is Codim< codim > :: EntityType.
          In general these types are the same, but for overloaded entities on grid parts
          this can differ.
        */
      template <class Entity>
      const Entity& convert( const Entity& entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().convert( entity ) );
        return asImp().convert( entity );
      }

      /** @brief Return true if the given entity is contained in this grid view
       * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
       *
       * \note If the input element e is not an element of the grid, then
       *       the result of contains() is undefined.
       */
      template<class EntityType>
      bool contains (const EntityType& e) const
      {
        return asImp().contains(e);
      }

      int overlapSize ( int codim ) const
      {
        DUNE_THROW( NotImplemented, "Method overlapSize() not implemented yet" );
      }

      int ghostSize( int codim ) const
      {
        DUNE_THROW( NotImplemented, "Method ghostSize() not implemented yet" );
      }

      /** \brief return sequence number to update structures depending on the grid part
       *  \note The default returns DofManager< Grid > :: sequence ()
       */
      [[deprecated("Use DofManager::sequence instead!")]]
      int sequence () const
      {
        DUNE_THROW(NotImplemented,"GridPart::sequence is not part of the interface, use DofManager::sequence instead!");
        return -1;
      }

    protected:
      //! do not create explicit instances of this class
      GridPartInterface () = default;

    private:
      GridPartType &asImp () { return static_cast< GridPartType & >( *this ); }
      const GridPartType &asImp () const { return static_cast< const GridPartType & >( *this ); }
    };



    //! \brief Default implementation for the GridPart classes
    template< class GridPartTraits >
    class GridPartDefault
    : public GridPartInterface< GridPartTraits >
    {
      typedef GridPartDefault< GridPartTraits > ThisType;
      typedef GridPartInterface< GridPartTraits > BaseType;

    public:
      //! \brief Type of the Traits
      typedef GridPartTraits Traits;
      //! \brief Grid implementation
      typedef typename Traits::GridType GridType;
      //! \brief Index set implementation
      typedef typename Traits::IndexSetType IndexSetType;

      //! type of intersection iterator
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

      //! type of intersection
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      //! \brief Collective communication
      typedef typename Traits::CommunicationType CommunicationType;

      //! \brief type of DofManager
      typedef DofManager< GridType >  DofManagerType;

      //! type of boundary id provider specialized for each grid type
      typedef BoundaryIdProvider< GridType > BoundaryIdProviderType;

    protected:
      GridType       *grid_;
      DofManagerType *dofManager_;

      //! constructor
      GridPartDefault ( GridType &grid )
      : grid_( &grid ),
        dofManager_( &DofManagerType :: instance( this->grid() ) )
      {}

      GridPartDefault ( const ThisType &other )
      : grid_( other.grid_ ),
        dofManager_( &DofManagerType :: instance( this->grid() ) )
      {}

      GridPartDefault& operator= ( const ThisType &other )
      {
        grid_ = other.grid_;
        dofManager_ = &DofManagerType :: instance( grid() );
        return *this;
      }

    public:
      //! Returns const reference to the underlying grid
      const GridType &grid () const { assert( grid_ ); return *grid_; }

      //! Returns reference to the underlying grid
      GridType &grid () { assert( grid_ ); return *grid_; }

      //! \brief returns true if current grid part represents a conforming grid
      bool isConforming() const
      {
        // default is to return the static variant
        return bool(BaseType::conforming);
      }

      /** \brief obtain number of entities in a given codimension */
      int size ( int codim ) const
      {
        return BaseType::indexSet().size( codim );
      }

      /** \brief obtain number of entities with a given geometry type */
      int size ( const GeometryType &type ) const
      {
        return BaseType::indexSet().size( type );
      }

      template<class EntityType>
      bool contains (const EntityType& e) const
      {
        return BaseType::indexSet().contains(e);
      }

      /** \brief obtain collective communication object */
      const CommunicationType &comm () const
      {
        return grid().comm();
      }

      /** \brief \copydoc GridPartInterface::entity
       *
       * \tparam  EntitySeed  entity seed from which to create entity
       *
       * The default implementation simply forwards to the corresponding
       * method on the grid.
       */
      template < class EntitySeed >
      typename Traits::template Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return grid().entity( seed );
      }

      /** \brief \copydoc GridPartInterface::convert

          \note  The default implementation does nothing but return the same entity
       */
      template <class Entity>
      const Entity& convert( const Entity& entity ) const
      {
        return entity;
      }

      //! \brief \copydoc GridPartInterface::entity
      [[deprecated("Use BoundaryIdProvider::boundaryId directly!")]]
      int boundaryId ( const IntersectionType &intersection ) const
      {
        return BoundaryIdProviderType::boundaryId( intersection );
      }
    };

  /** @} */

    template< class Entity >
    struct GridEntityAccess;

    template< int codim, int dim, class Grid, template< int, int, class > class EntityImpl >
    struct GridEntityAccess< Dune::Entity< codim, dim, Grid, EntityImpl > >
    {
      typedef Dune::Entity< codim, dim, Grid, EntityImpl > EntityType;
      typedef Dune::Entity< codim, dim, Grid, EntityImpl > GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return entity;
      }
    };

    template< class Entity >
    const typename GridEntityAccess< Entity >::GridEntityType &
    gridEntity ( const Entity &entity )
    {
      return GridEntityAccess< Entity >::gridEntity( entity );
    }

  } // namespace Fem

} // namespace Dune

#endif // #define DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH
