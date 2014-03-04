#ifndef DUNE_FEM_GRIDPART_LEAFGRIDPART_HH
#define DUNE_FEM_GRIDPART_LEAFGRIDPART_HH

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
    struct LeafGridPartTraits;


    // LeafGridPart
    // ------------

    //! \brief Selects the leaf level of a grid
    template< class GridImp >
    class LeafGridPart
    : public GridPartDefault< LeafGridPartTraits< GridImp > >
    {
      typedef LeafGridPart< GridImp > ThisType;
      typedef GridPartDefault< LeafGridPartTraits< GridImp > > BaseType;

    public:
      //- Public typedefs and enums
      //! Type definitions
      typedef LeafGridPartTraits< GridImp > Traits;

      //! Grid implementation type
      typedef typename Traits::GridType GridType;
      //! The leaf index set of the grid implementation
      typedef typename Traits::IndexSetType IndexSetType;
      
      //! The corresponding IntersectionIterator 
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

      typedef typename IntersectionIteratorType::Intersection IntersectionType;
    
      //! the leaf grid view from the grid 
      typedef typename GridType::template Partition< All_Partition >::LeafGridView LeafGridView;

    private:
      typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

    public:
      //- Public methods
      //! Constructor
      explicit LeafGridPart ( GridType &grid )
      : BaseType( grid ),
        leafGridView_( grid.leafGridView() ),
        isetWrapper_( grid )
      {}

      //! copy constructor
      LeafGridPart ( const ThisType &other )
      : BaseType( other ),
        leafGridView_( other.leafGridView_ ),
        isetWrapper_( other.grid() )
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
        return leafGridView_.template begin< codim, pitype >();
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename BaseType::template Codim< codim >::IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return leafGridView_.template end< codim, pitype >();
      }

      //! ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType ibegin ( const EntityCodim0Type &entity ) const
      {
        return leafGridView_.ibegin( entity );
      }
      
      //! iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend ( const EntityCodim0Type &entity ) const
      {
        return leafGridView_.iend( entity );
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return intersection.boundaryId();
      }

      //! Returns maxlevel of the grid
      int level() const { return grid().maxLevel(); }

      //! corresponding communication method for this grid part
      template <class DataHandleImp,class DataType>
      void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                       InterfaceType iftype, CommunicationDirection dir) const 
      {
        leafGridView_.communicate( data, iftype, dir );
      }

    private: 
      //! leaf grid view 
      LeafGridView leafGridView_ ;
      //! GridDefaultIndexSet Wrapper 
      IndexSetType isetWrapper_;
    };



    // LeafGridPartTraits
    // ------------------

    //! Type definitions for the LeafGridPart class
    template< class GridImp >
    struct LeafGridPartTraits
    {
      /** \brief The type of the grid */
      typedef GridImp GridType;

      /** \brief The type of the corresponding grid part class */
      typedef LeafGridPart< GridImp > GridPartType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;

      /** \brief The appropriate index set */
      typedef WrappedLeafIndexSet<GridType> IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      static const InterfaceType indexSetInterfaceType = All_All_Interface;

      /** \brief The appropriate intersection iterator */
      typedef typename GridType::template Codim<0>::Entity::LeafIntersectionIterator
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
          typedef typename GridType::template Codim< codim >::template Partition< pitype >::LeafIterator
            IteratorType;
        };
      };

      typedef typename GridType::CollectiveCommunication CollectiveCommunicationType;

      //! \brief is true if grid on this view only has conforming intersections 
      static const bool conforming = Dune::Capabilities::isLeafwiseConforming< GridType >::v;
    };



    // Leaf grid part capabilities
    // ---------------------------

    namespace GridPartCapabilities
    {

      template< class GridType >
      struct hasGrid< LeafGridPart< GridType > >
      {
        static const bool v = true;
      };

      template< class GridType >
      struct hasSingleGeometryType< LeafGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId;
      };

      template< class GridType >
      struct isCartesian< LeafGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isCartesian< GridType >::v;
      };

      template< class GridType, int codim  >
      struct hasEntity< LeafGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< GridType, codim >::v;
      };

      template< class GridType >
      struct isParallel< LeafGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isParallel< GridType >::v;
      };

      template< class GridType, int codim >
      struct canCommunicate< LeafGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< GridType, codim >::v;
      };

      template< class GridType >
      struct isConforming< LeafGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< GridType >::v;
      };

    } // end namespace GridPartCapabilities

    /** @} */

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: LeafGridPart ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_LEAFGRIDPART_HH
