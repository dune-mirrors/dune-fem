#ifndef DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH
#define DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH

#include <utility>
#include <functional>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridView, class Implementation, bool storeCopy=true >
    class GridView2GridPart;



#ifndef DOXYGEN

    // GridView2GridPartTraits
    // -----------------------

    template< class GridView, class Implementation, bool storeCopy >
    struct GridView2GridPartTraits
    {
      typedef Implementation GridPartType;

      typedef GridView GridViewType;
      static const bool conforming = GridView::conforming;

      typedef typename GridViewType::Grid GridType;
      typedef typename GridViewType::Communication CommunicationType;

      typedef typename GridView::IndexSet IndexSetType;

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

    template< class GridView, class Implementation, bool storeCopy >
    class GridView2GridPart
      : public GridPartInterface< GridView2GridPartTraits< GridView, Implementation, storeCopy > >
    {
      typedef GridView2GridPart< GridView, Implementation, storeCopy > ThisType;
      typedef GridView2GridPartTraits< GridView, Implementation, storeCopy > TraitsType;
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

      /** \copydoc Dune::Fem::GridPartInterface::CommunicationType */
      typedef typename BaseType::CommunicationType CommunicationType;

    private:
      typedef DofManager< GridType > DofManagerType;

      auto initGv( const GridView &gridView )
      {
        if constexpr ( storeCopy )
          return gridView;
        else
          return &gridView;
      }
    public:
      using BaseType::grid;
      using BaseType::boundaryId;

      /** \name Construction
       *  \{
       */

      explicit GridView2GridPart ( const GridView &gridView )
        : gridView_( initGv( gridView ) ),
          indexSet_( &this->gridView().indexSet() )
      {}

      explicit GridView2GridPart ( GridView &&gridView )
        : gridView_( std::move( gridView ) ),
          indexSet_( &this->gridView().indexSet() )
      {
        // this should not be called if we only store a pointer
        assert( storeCopy );
      }

      GridView2GridPart ( const ThisType &rhs )
        : gridView_( rhs.gridView_ ),
          indexSet_( &gridView().indexSet() )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::GridPartInterface::grid */
      const GridType &grid () const { return gridView().grid(); }
      GridType &grid () { return const_cast< GridType & >( gridView().grid() ); }

      /** \copydoc Dune::Fem::GridPartInterface::indexSet */
      const IndexSetType &indexSet () const { assert( indexSet_ ); return *indexSet_; }

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
        return gridView().template begin< codim, pitype >();
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
        return gridView().template end< codim, pitype >();
      }

      /** \copydoc Dune::Fem::GridPartInterface::ibegin */
      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return gridView().ibegin( entity );
      }

      /** \copydoc Dune::Fem::GridPartInterface::iend */
      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return gridView().iend( entity );
      }

      /** \copydoc Dune::Fem::GridPartInterface::comm */
      const CommunicationType &comm () const { return gridView().comm(); }

      /** \copydoc Dune::Fem::GridPartInterface::communicate */
      template< class DataHandle, class DataType >
      void communicate ( CommDataHandleIF< DataHandle, DataType > &dataHandle,
                         InterfaceType interface, CommunicationDirection direction ) const
      {
        gridView().communicate( dataHandle, interface, direction );
      }

      /** \copydoc Dune::Fem::GridPartInterface::sequence */
      [[deprecated("Use DofManager::sequence instead!")]]
      int sequence () const { return -1; }

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
      operator const GridView& () const { return gridView(); }

      /** \brief return reference to internal grid view */
      const GridView& gridView() const
      {
        if constexpr ( storeCopy )
        {
          return gridView_;
        }
        else
        {
          assert( gridView_ );
          return *gridView_;
        }
      }

      /** \} */

    protected:
      template< int codim >
      const typename Codim< codim >::EntityType &
      convert( const typename Codim< codim >::EntityType &entity ) const
      {
        return entity;
      }

      std::conditional_t<storeCopy, GridView, const GridView* > gridView_;
      const IndexSetType* indexSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_GRIDVIEW2GRIDPART_HH
