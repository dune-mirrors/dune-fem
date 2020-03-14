#ifndef DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH
#define DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH

#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/nonadaptiveindexset.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridView, class Implementation >
    class GridView2GridPart;



#ifndef DOXYGEN

    // GridView2GridPartTraits
    // -----------------------

    template< class GridView, class Implementation >
    struct GridView2GridPartTraits
    {
      typedef Implementation GridPartType;

      typedef GridView GridViewType;
      static const bool conforming = GridView::conforming;

      typedef typename GridViewType::Grid GridType;
      typedef typename GridViewType::CollectiveCommunication CollectiveCommunicationType;

      typedef NonAdaptiveIndexSet< typename GridView::IndexSet > IndexSetType;

      template< int codim >
      struct Codim
      {
        typedef typename GridViewType::template Codim< codim >::Entity EntityType;
        typedef typename GridType::template Codim< codim >::EntitySeed EntitySeedType;

        typedef typename GridViewType::template Codim< codim >::Geometry GeometryType;
        typedef typename GridViewType::template Codim< codim >::LocalGeometry LocalGeometryType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridViewType::template Codim< codim >::template Partition< pitype >::Iterator IteratorType;
        };
      };

      typedef typename GridViewType::IntersectionIterator IntersectionIteratorType;

      typedef TwistUtility< GridType > TwistUtilityType;

      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      static const InterfaceType indexSetInterfaceType = All_All_Interface;
    };

#endif // #ifndef DOXYGEN



    // GridView2GridPart
    // -----------------

    template< class GridView, class Implementation >
    class GridView2GridPart
      : public GridPartInterface< GridView2GridPartTraits< GridView, Implementation > >
    {
      typedef GridView2GridPart< GridView, Implementation > ThisType;
      typedef GridView2GridPartTraits< GridView, Implementation > TraitsType;
      typedef GridPartInterface< TraitsType > BaseType;

    public:
      /** \copydoc Dune::Fem::GridPartInterface::GridType */
      typedef typename BaseType::GridType GridType;

      /** \copydoc Dune::Fem::GridPartInterface::GridViewType */
      typedef typename BaseType::GridViewType GridViewType;

      template< int codim >
      struct Codim
        : public BaseType::template Codim< codim >
      {};

      /** \copydoc Dune::Fem::GridPartInterface::IntersectionIteratorType */
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;

      /** \copydoc Dune::Fem::GridPartInterface::IndexSetType */
      typedef typename BaseType::IndexSetType IndexSetType;

      /** \copydoc Dune::Fem::GridPartInterface::CollectiveCommunicationType */
      typedef typename BaseType::CollectiveCommunicationType CollectiveCommunicationType;

    private:
      typedef DofManager< GridType > DofManagerType;

    public:
      using BaseType::grid;
      using BaseType::boundaryId;

      /** \name Construction
       *  \{
       */

      explicit GridView2GridPart ( const GridView &gridView )
        : gridView_( gridView ),
          indexSet_( gridView_.indexSet() ),
          dofManager_( DofManagerType::instance( gridView_.grid() ) )
      {}

      explicit GridView2GridPart ( GridView &&gridView )
        : gridView_( std::move( gridView ) ),
          indexSet_( gridView_.indexSet() ),
          dofManager_( DofManagerType::instance( gridView_.grid() ) )
      {}

      GridView2GridPart ( const ThisType &rhs )
        : gridView_( rhs.gridView_ ),
          indexSet_( gridView_.indexSet() ),
          dofManager_( DofManagerType::instance( rhs.grid() ) )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      const GridType &grid () const { return gridView_.grid(); }

      /** \copydoc Dune::Fem::GridPartInterface::indexSet */
      const IndexSetType &indexSet () const { return indexSet_; }

      /** \copydoc Dune::Fem::GridPartInterface::begin */
      template< int codim >
      typename Codim< codim >::IteratorType begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      /** \copydoc Dune::Fem::GridPartInterface::begin */
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType begin () const
      {
        return gridView_.template begin< codim, pitype >();
      }

      /** \copydoc Dune::Fem::GridPartInterface::end */
      template< int codim >
      typename Codim< codim >::IteratorType end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      /** \copydoc Dune::Fem::GridPartInterface::end */
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType end () const
      {
        return gridView_.template end< codim, pitype >();
      }

      /** \copydoc Dune::Fem::GridPartInterface::ibegin */
      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return gridView_.ibegin( entity );
      }

      /** \copydoc Dune::Fem::GridPartInterface::iend */
      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return gridView_.iend( entity );
      }

      /** \copydoc Dune::Fem::GridPartInterface::comm */
      const CollectiveCommunicationType &comm () const { return gridView_.comm(); }

      /** \copydoc Dune::Fem::GridPartInterface::communicate */
      template< class DataHandle, class DataType >
      void communicate ( CommDataHandleIF< DataHandle, DataType > &dataHandle,
                         InterfaceType interface, CommunicationDirection direction ) const
      {
        gridView_.communicate( dataHandle, interface, direction );
      }

      /** \copydoc Dune::Fem::GridPartInterface::sequence */
      int sequence () const { return dofManager_.sequence(); }

      /** \copydoc Dune::Fem::GridPartInterface::entity */
      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return grid().entity( seed );
      }

      /** \copydoc Dune::Fem::GridPartInterface::sequence */
      template <class Entity>
      const Entity &convert( const Entity& entity ) const
      {
        return convert< Entity::codimension >( entity );
      }

      /** \brief cast to underlying grid view */
      explicit operator GridView () const { return gridView_; }

      /** \} */

    private:
      template< int codim >
      const typename Codim< codim >::EntityType &
      convert( const typename Codim< codim >::EntityType &entity ) const
      {
        return entity;
      }

      GridView gridView_;
      IndexSetType indexSet_;
      DofManagerType &dofManager_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH
