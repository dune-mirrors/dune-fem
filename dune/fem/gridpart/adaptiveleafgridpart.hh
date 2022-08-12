#ifndef DUNE_FEM_ADAPTIVELEAFGRIDPART_HH
#define DUNE_FEM_ADAPTIVELEAFGRIDPART_HH

//- dune-grid includes
#include <dune/grid/common/gridview.hh>

//- dune-fem includes
#include <dune/fem/gridpart/adaptiveleafindexset.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/storage/singletonlist.hh>

namespace Dune
{

  namespace Fem
  {

    /*- see dune/grid/common/gridenums.hh
      enum InterfaceType {
            InteriorBorder_InteriorBorder_Interface=0, //!< send/receive interior and border entities
            InteriorBorder_All_Interface=1,            //!< send interior and border, receive all entities
            Overlap_OverlapFront_Interface=2,          //!< send overlap, receive overlap and front entities
            Overlap_All_Interface=3,                   //!< send overlap, receive all entities
            All_All_Interface=4                        //!< send all and receive all entities
      };
      enum PartitionIteratorType {
            Interior_Partition=0,       //!< only interior entities
            InteriorBorder_Partition=1, //!< interior and border entities
            Overlap_Partition=2,        //!< only overlap entities
            OverlapFront_Partition=3,   //!< overlap and front entities
            All_Partition=4,            //!< all entities
            Ghost_Partition=5           //!< only ghost entities
      };
    */

    template <PartitionIteratorType ittype>
    struct IteratorToInterface
    {
      static const InterfaceType value = InteriorBorder_All_Interface;
    };
    template <>
    struct IteratorToInterface< InteriorBorder_Partition >
    {
      static const InterfaceType value = InteriorBorder_InteriorBorder_Interface;
    };

    /////////////////////////////////////////////////////////////////////////
    //
    //  --AdaptiveLeafIndexGridPart
    //
    /////////////////////////////////////////////////////////////////////////

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations.

        The underlying \ref AdaptiveLeafIndexSet "index set" is defined for
        entities of all codimensions.
    */
    template< class TraitsImp >
    class AdaptiveGridPartBase
    : public GridPartDefault< TraitsImp >
    {
      typedef AdaptiveGridPartBase< TraitsImp > ThisType;
      typedef GridPartDefault< TraitsImp > BaseType;

    public:
      //! Type definitions
      typedef TraitsImp Traits;

      //! Grid implementation type
      typedef typename Traits :: GridPartType GridPartType;
      //! Grid implementation type
      typedef typename Traits :: GridType GridType;
      //! The leaf index set of the grid implementation
      typedef typename Traits :: IndexSetType IndexSetType;

      //! Struct providing types of the leaf iterators on codimension codim
      template< int codim >
      struct Codim
      : public BaseType :: template Codim< codim >
      {};

    private:
      typedef typename GridType::LeafGridView LeafGridView;

    public:
      //! type of intersection iterator
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

      //! type of intersection
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      typedef std::integral_constant< bool, false > NoIndexSetType;

      typedef GridPartType GridViewType;

    protected:
      // key type for singleton list is grid pointer
      typedef SingletonList < const GridType*, IndexSetType > IndexSetProviderType;

      // type of entity with codimension zero
      typedef typename Codim< 0 > :: EntityType ElementType;

      // the leaf grid view
      LeafGridView leafGridView_ ;

      // reference to index set
      std::shared_ptr< IndexSetType > indexSet_;

    public:
      //! constructor
      explicit AdaptiveGridPartBase ( GridType &grid )
      : BaseType( grid ),
        leafGridView_( grid.leafGridView() ),
        indexSet_( &IndexSetProviderType::getObject( &grid ),
                   typename IndexSetProviderType::Deleter() )
      {}

      //! Copy Constructor
      AdaptiveGridPartBase ( const ThisType &other ) = default;

      AdaptiveGridPartBase& operator= ( const AdaptiveGridPartBase& other ) = default;

    protected:
      //! Constructor constructing object held by index set (for iterator access)
      AdaptiveGridPartBase ( GridType& grid, const NoIndexSetType& noIndexSet )
      : BaseType( grid ),
        leafGridView_( grid.leafGridView() ),
        indexSet_() // not created because noIndexSet was passed
      {}

    public:
      using BaseType::grid;

      //! Returns reference to index set of the underlying grid
      const IndexSetType &indexSet () const
      {
        assert( indexSet_ );
        return *indexSet_;
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim > :: IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      //! Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: IteratorType
      begin () const
      {
        return leafGridView_.template begin< codim, pitype >();
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim > :: IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: IteratorType
      end () const
      {
        return leafGridView_.template end< codim, pitype >();
      }

      //! ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType
      ibegin ( const ElementType &entity ) const
      {
        return leafGridView_.ibegin( entity );
      }

      //! iend of corresponding intersection iterator for given entity
      IntersectionIteratorType
      iend ( const ElementType &entity ) const
      {
        return leafGridView_.iend( entity );
      }

      //! corresponding communication method for this grid part
      template< class DataHandle, class Data >
      decltype( auto ) communicate ( CommDataHandleIF< DataHandle, Data > &data, InterfaceType iftype, CommunicationDirection dir ) const
      {
        return leafGridView_.communicate( data, iftype, dir );
      }

    protected:
      const GridPartType& asImp() const
      {
        return static_cast<const GridPartType &> (*this);
      }

      GridPartType& asImp()
      {
        return static_cast<GridPartType &> (*this);
      }
    };

    /** @addtogroup AdaptiveLeafGP
        \brief GridPart for Dune::AdaptiveLeafIndexSet.
        The underlying index set is
        singleton for each grid object.
        Uses very efficient index sets specially
        designed for problems with constantly changing underlying grid.
    */
    template< class Grid, PartitionIteratorType idxpitype = All_Partition, bool onlyCodimensionZero = false >
    class AdaptiveLeafGridPart;

    //! Type definitions for the LeafGridPart class
    template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero >
    class AdaptiveLeafGridPartTraits
    {
    public:
      //! type of the grid
      typedef Grid GridType;

      //! type of the grid part , i.e. this type
      typedef AdaptiveLeafGridPart< GridType, idxpitype, onlyCodimensionZero > GridPartType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;

      typedef typename GridType::Communication CommunicationType;

    protected:
      // choose the AdaptiveIndexSet (based on the HierarchicIndexSet)
      // to be revised
      template < int dummy, bool onlyCodimZero >
      struct AdaptiveLeafIndexSetChooser
      {
        static const PartitionIteratorType indexSetPartitionType = idxpitype;
        static const InterfaceType indexSetInterfaceType = IteratorToInterface<idxpitype>::value;
        typedef AdaptiveLeafIndexSet< GridPartType > IndexSetType;
      };

      template <int dummy>
      struct AdaptiveLeafIndexSetChooser<dummy, true >
      {
        static const PartitionIteratorType indexSetPartitionType = idxpitype;
        static const InterfaceType indexSetInterfaceType = IteratorToInterface<idxpitype>::value;
        typedef DGAdaptiveLeafIndexSet< GridPartType > IndexSetType;
      };

      // also for Cartesian grids (e.g. YaspGrid) use adaptive leaf index set in parallel
      typedef AdaptiveLeafIndexSetChooser<-1, onlyCodimensionZero > IndexSetChooserType;

    public:
      //! type of the index set
      typedef typename IndexSetChooserType::IndexSetType IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = IndexSetChooserType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = IndexSetChooserType::indexSetInterfaceType;

      // type of intersection iterator
      typedef typename GridType::LeafGridView::IntersectionIterator IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridType::template Codim< codim >::Geometry      GeometryType;
        typedef typename GridType::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridType::template Codim< codim >::Entity        EntityType;
        typedef typename GridType::template Codim< codim >::EntitySeed    EntitySeedType;

        // GridView typedefs interface
        typedef GeometryType       Geometry;
        typedef LocalGeometryType  LocalGeometry;
        typedef EntityType         Entity;
        typedef EntitySeedType     EntitySeed;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridType::template Codim< codim >::template Partition< pitype >::LeafIterator IteratorType;
          // GridView typedef
          typedef IteratorType Iterator;
        };

        typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;
        // GridView typedef
        typedef IteratorType Iterator;
      };

      //! \brief is true if grid on this view only has conforming intersections
      static const bool conforming = Dune::Capabilities::isLeafwiseConforming< GridType > :: v;
    };

    template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero >
    class AdaptiveLeafGridPart
      : public AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > >
    {
      typedef AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > > BaseType;
    public:
      typedef typename BaseType :: NoIndexSetType  NoIndexSetType;
      typedef typename BaseType :: GridType        GridType;
      typedef typename BaseType :: GridViewType    GridViewType;
      typedef typename BaseType :: GridPartType    GridPartType;
      //! Constructor
      explicit AdaptiveLeafGridPart ( GridType &grid )
      : BaseType( grid )
      {}

      //! copy constructor (for construction from IndexSet, no public use)
      AdaptiveLeafGridPart ( GridType& grid, const NoIndexSetType& dummy )
      : BaseType( grid, dummy )
      {}

      //! copy constructor
      AdaptiveLeafGridPart ( const AdaptiveLeafGridPart& other ) = default;
    };

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations.

        The underlying \ref DGAdaptiveLeafIndexSet "index set" is defined
        only for codimension 0.
    */
    template< class Grid, PartitionIteratorType idxpitype = All_Partition >
    using DGAdaptiveLeafGridPart = AdaptiveLeafGridPart< Grid, idxpitype, true >;


    template< class Grid, PartitionIteratorType idxpitype = All_Partition >
    class IntersectionAdaptiveLeafGridPart ;

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations including indices for intersections.

        The underlying \ref IntersectionAdaptiveLeafIndexSet "index set" is defined
        also for codimension -1 (intersections).
    */
    template< class Grid, PartitionIteratorType idxpitype >
    class IntersectionAdaptiveLeafGridPartTraits : public AdaptiveLeafGridPartTraits< Grid, idxpitype, false>
    {
    public:
      //! type of the grid part , i.e. this type
      typedef IntersectionAdaptiveLeafGridPart< Grid, idxpitype > GridPartType;

      //! type of the index set
      typedef IntersectionAdaptiveLeafIndexSet< GridPartType > IndexSetType;
    };

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations including indices for intersections.

        The underlying \ref IntersectionAdaptiveLeafIndexSet "index set" is defined
        also for codimension -1 (intersections).
    */
    template< class Grid, PartitionIteratorType idxpitype >
    class IntersectionAdaptiveLeafGridPart
    : public AdaptiveGridPartBase< IntersectionAdaptiveLeafGridPartTraits< Grid, idxpitype > >
    {
      typedef AdaptiveGridPartBase< IntersectionAdaptiveLeafGridPartTraits< Grid, idxpitype > > BaseType;
    public:
      typedef typename BaseType :: NoIndexSetType  NoIndexSetType;
      typedef typename BaseType :: GridType GridType;
      //! Constructor
      explicit IntersectionAdaptiveLeafGridPart( GridType &grid )
      : BaseType( grid )
      {
      }

      //! copy constructor (for construction from IndexSet, no public use)
      IntersectionAdaptiveLeafGridPart ( GridType& grid, const NoIndexSetType& noIndexSet )
      : BaseType( grid, noIndexSet )
      {
      }

      //! copy constructor
      IntersectionAdaptiveLeafGridPart( const IntersectionAdaptiveLeafGridPart& other ) = default;
    };



    // Capabilities
    // ------------

    namespace GridPartCapabilities
    {

      // Capabilities for AdaptiveLeafGridPart
      // -------------------------------------

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct hasGrid< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct hasSingleGeometryType< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct isCartesian< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero, int codim  >
      struct hasEntity< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero, int codim >
      struct canCommunicate< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct isConforming< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };


      // Capabilities for DGAdaptiveLeafGridPart
      // ---------------------------------------

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasGrid< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasSingleGeometryType< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isCartesian< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim  >
      struct hasEntity< DGAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim >
      struct canCommunicate< DGAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isConforming< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };


      // Capbilities for IntersectionAdaptiveLeafGridPart
      // ------------------------------------------------

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasGrid< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasSingleGeometryType< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isCartesian< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim  >
      struct hasEntity< IntersectionAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim >
      struct canCommunicate< IntersectionAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isConforming< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ADAPTIVELEAFGRIDPART_HH
