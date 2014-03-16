#ifndef DUNE_FEM_GRIDPART_LEVELGRIDPART_HH
#define DUNE_FEM_GRIDPART_LEVELGRIDPART_HH

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>


namespace Dune
{

  namespace Fem
  {

    /** 
     * @addtogroup GridPart
     *
     * @{ 
     */

    // Forward declarations
    // --------------------

    template< class >
    class GridPartDefault;
    template< class >
    class LevelGridPartTraits;



    // LevelGridPart
    // -------------

    //! \brief Selects a specific level of a grid
    template< class GridImp >
    class LevelGridPart
    : public GridPartDefault< LevelGridPartTraits< GridImp > >
    {
      typedef LevelGridPart< GridImp > ThisType;
      typedef GridPartDefault< LevelGridPartTraits< GridImp > > BaseType;

    public:
      //- Public typedefs and enums
      //! Corresponding type definitions
      typedef LevelGridPartTraits< GridImp > Traits;

      //! Grid implementation
      typedef typename Traits::GridType GridType;
      //! Level index set that corresponds to the grid
      typedef typename Traits::IndexSetType IndexSetType;
      
      //! The corresponding IntersectionIterator 
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      typedef typename GridType::template Partition< All_Partition > :: LevelGridView  LevelGridView;

    private:
      typedef typename GridType::template Codim< 0 >::Entity EntityCodim0Type;

    public:
      //! Constructor
      LevelGridPart ( GridType &grid, const int level )
      : BaseType( grid ),
        levelGridView_( grid.levelGridView( level ) ),
        isetWrapper_( grid, level ),
        level_( level )
      {}
      
      //! Constructor, choosing maxLevel
      explicit LevelGridPart ( GridType &grid )
      : BaseType( grid ),
        levelGridView_( grid.levelGridView( grid.maxLevel() ) ),
        isetWrapper_( grid, grid.maxLevel() ),
        level_( grid.maxLevel() )
      {}
      
      //! copy constructor
      LevelGridPart ( const ThisType &other )
      : BaseType( other ),
        levelGridView_( other.levelGridView_ ),
        isetWrapper_( other.grid(), other.level_ ),
        level_( other.level_ )
      {}

      using BaseType::grid;

      //! Returns reference to index set of the underlying grid
      const IndexSetType &indexSet () const
      {
        return isetWrapper_;
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename BaseType::template Codim< codim >::IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      //! Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        return levelGridView_.template begin< codim, pitype >();
      }

      //! Begin iterator on the GridPart's level
      template< int codim >
      typename BaseType::template Codim< codim >::IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! End iterator on the GridPart's level
      template< int codim, PartitionIteratorType pitype >
      typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return levelGridView_.template end< codim, pitype >();
      }

      //! ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType ibegin ( const EntityCodim0Type &entity ) const
      {
        return levelGridView_.ibegin( entity );
      }
      
      //! iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend ( const EntityCodim0Type &entity ) const
      {
        return levelGridView_.iend( entity );
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return intersection.boundaryId();
      }

      //! Level which this GridPart belongs to
      int level() const { return level_; }

      //! corresponding communication method for this grid part
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        levelGridView_.communicate( data, iftype, dir );
      }

    private:
      LevelGridView levelGridView_;
      //! GridDefaultIndexSet Wrapper 
      IndexSetType isetWrapper_;
      const int level_;
    };



    // LevelGridPartTraits
    // -------------------

    //! Type definitions for the LevelGridPart class
    template< class GridImp >
    struct LevelGridPartTraits
    {
      /** \brief The type of the grid */
      typedef GridImp GridType;

      /** \brief The type of the corresponding grid part class */
      typedef LevelGridPart< GridImp > GridPartType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;

      /** \brief The appropriate index set */
      typedef WrappedLevelIndexSet<GridType> IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      static const InterfaceType indexSetInterfaceType = All_All_Interface;

        /** \brief The appropriate intersection iterator */
      typedef typename GridType::template Codim<0>::Entity::LevelIntersectionIterator
        IntersectionIteratorType;

        /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      template< int codim >
      struct Codim
      {
        typedef typename GridType::template Codim< codim >::Geometry GeometryType;
        typedef typename GridType::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridType::template Codim< codim >::Entity EntityType;
        typedef typename GridType::template Codim< codim >::EntityPointer EntityPointerType;

        typedef typename GridType::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridType::template Codim< codim >::template Partition< pitype >::LevelIterator IteratorType;
        };
      };

      typedef typename GridType::CollectiveCommunication CollectiveCommunicationType;

      //! \brief is true if grid on this view only has conforming intersections
      static const bool conforming = Dune::Capabilities::isLevelwiseConforming< GridType >::v;
    };

    /** @} */

    // Level grid part capabilities
    // ----------------------------

    namespace GridPartCapabilities
    {

      template< class GridType >
      struct hasGrid< LevelGridPart< GridType > >
      {
        static const bool v = true;
      };

      template< class GridType >
      struct hasSingleGeometryType< LevelGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId;
      };

      template< class GridType >
      struct isCartesian< LevelGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isCartesian< GridType >::v;
      };

      template< class GridType, int codim  >
      struct hasEntity< LevelGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< GridType, codim >::v;
      };

      template< class GridType >
      struct isParallel< LevelGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isParallel< GridType >::v;
      };

      template< class GridType, int codim >
      struct canCommunicate< LevelGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< GridType, codim >::v;
      };

      template< class GridType >
      struct isConforming< LevelGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isLevelwiseConforming< GridType >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: LevelGridPart ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_LEVELGRIDPART_HH
