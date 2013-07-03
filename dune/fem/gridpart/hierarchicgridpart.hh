#ifndef DUNE_FEM_GRIDPART_HIERARCHICGRIDPART_HH
#define DUNE_FEM_GRIDPART_HIERARCHICGRIDPART_HH

//- dune-fem includes
#include <dune/fem/version.hh>
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
    class HierarchicGridPartTraits;



    // HierarchicGridPart
    // ------------------

    /** \brief Selects the leaf level of a grid together with the 
        HierarchicIndexSet available for ALUGrid and AlbertaGrid. 
        The HierarchicIndexSet is basically the LocalIdSet of the grid 
        extended by a size method to implement the IndexSet interface. 
        For all other grids the default LeafIndexSet is selected.
    */
    template< class GridImp >
    class HierarchicGridPart
    : public GridPartDefault< HierarchicGridPartTraits< GridImp > >
    {
      typedef HierarchicGridPart< GridImp > ThisType;
      typedef GridPartDefault< HierarchicGridPartTraits< GridImp > > BaseType;

    public:
      //- Public typedefs and enums
      //! Type definitions
      typedef HierarchicGridPartTraits< GridImp > Traits;

      //! Grid implementation type
      typedef typename Traits::GridType GridType;
      //! The leaf index set of the grid implementation
      typedef typename Traits::IndexSetType IndexSetType;
      
      //! The corresponding Intersection
      typedef typename Traits::IntersectionType IntersectionType ;

      //! The corresponding IntersectionIterator 
      typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

      //! the leaf grid view from the grid 
      typedef typename GridType :: template Partition < All_Partition > :: LeafGridView LeafGridView;

    private:
      typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

    public:
      //- Public methods
      //! constructor 
      DUNE_VERSION_DEPRECATED(1,4,remove)
      explicit HierarchicGridPart ( GridType &grid )
      : BaseType( grid ),
        leafView_( grid.leafView() ),
        isetWrapper_( grid )
      {}

      //! constructor
      DUNE_VERSION_DEPRECATED(1,4,remove)
      HierarchicGridPart ( GridType &grid, const IndexSetType & )
      : BaseType( grid ),
        leafView_( grid.leafView() ),
        isetWrapper_( grid )
      {}

      //! copy constructor
      DUNE_VERSION_DEPRECATED(1,4,remove)
      HierarchicGridPart ( const ThisType &other )
      : BaseType( other ),
        leafView_( other.leafView_ ),
        isetWrapper_( other.grid() )
      {}

      //! Returns reference to index set of the underlying grid
      const IndexSetType &indexSet () const
      {
        return isetWrapper_;
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename BaseType :: template Codim< codim > :: IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      //! Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Traits :: template Codim< codim > :: template Partition< pitype > :: IteratorType
      begin () const
      {
        return leafView_.template begin< codim, pitype >();
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename BaseType :: template Codim< codim > :: IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Traits :: template Codim< codim > :: template Partition< pitype > :: IteratorType
      end () const
      {
        return leafView_.template end< codim, pitype >();
      }

      //! ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType ibegin(const EntityCodim0Type & en) const 
      {
        return en.ileafbegin();
      }
      
      //! iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend(const EntityCodim0Type & en) const 
      {
        return en.ileafend();
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return intersection.boundaryId();
      }

      //! Returns maxlevel of the grid
      int level() const { return this->grid().maxLevel(); }

      //! corresponding communication method for this grid part
      template <class DataHandleImp,class DataType>
      void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                       InterfaceType iftype, CommunicationDirection dir) const 
      {
        this->grid().communicate(data,iftype,dir);
      }

    private: 
      //! leaf grid view 
      LeafGridView leafView_ ;
      //! GridDefaultIndexSet Wrapper 
      IndexSetType isetWrapper_;
    };



    // HierarchicGridPartTraits
    // ------------------------

    //! Type definitions for the HierarchicGridPart class
    template< class GridImp >
    struct HierarchicGridPartTraits
    {
      /** \brief The type of the grid */
      typedef GridImp GridType;
      /** \brief The type of the corresponding grid part class */
      typedef HierarchicGridPart< GridImp > GridPartType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;

      /** \brief The appropriate index set */
      typedef WrappedHierarchicIndexSet<GridType> IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      static const InterfaceType indexSetInterfaceType = All_All_Interface;

        /** \brief The appropriate intersection */
      typedef typename GridType::Traits::
        LeafIntersection IntersectionType;

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
          typedef typename GridType::template Codim< codim >::template Partition< pitype >::LeafIterator IteratorType;
        };
      };

      //! \brief is true if grid on this view only has conforming intersections 
      static const bool conforming = Dune::Capabilities::isLeafwiseConforming< GridType >::v;
    };

    /** @} */



    // Hierarchic grid part capabilities
    // ---------------------------------

    namespace GridPartCapabilities
    {

      template< class GridType >
      struct hasGrid< HierarchicGridPart< GridType > >
      {
        static const bool v = true;
      };

      template< class GridType >
      struct hasSingleGeometryType< HierarchicGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId;
      };

      template< class GridType >
      struct isCartesian< HierarchicGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isCartesian< GridType >::v;
      };

      template< class GridType, int codim  >
      struct hasEntity< HierarchicGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< GridType, codim >::v;
      };

      template< class GridType >
      struct isParallel< HierarchicGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isParallel< GridType >::v;
      };

      template< class GridType, int codim >
      struct canCommunicate< HierarchicGridPart< GridType >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< GridType, codim >::v;
      };

      template< class GridType >
      struct isConforming< HierarchicGridPart< GridType > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< GridType >::v;
      };

    } // end namespace GridPartCapabilities

  } // end namespace Fem 

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: HierarchicGridPart ;
#endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_HIERARCHICGRIDPART_HH
