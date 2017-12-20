#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH

//- system includes
#include <cassert>
#include <memory>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridview.hh>

//- dune-fem includes
#include <dune/fem/gridpart/adaptiveleafindexset.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/fem/gridpart/filteredgridpart/capabilities.hh>
#include <dune/fem/gridpart/filteredgridpart/datahandle.hh>
#include <dune/fem/gridpart/filteredgridpart/intersection.hh>
#include <dune/fem/gridpart/filteredgridpart/intersectioniterator.hh>
#include <dune/fem/gridpart/filteredgridpart/iterator.hh>


namespace Dune
{

  namespace Fem
  {

    // Forward declarations
    // --------------------

    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet = false >
    class FilteredGridPart;



    // FilteredGridPartIndexSetSelector
    // --------------------------------

    template < class FilteredGP, class HostGP, bool useFilteredIndexSet >
    struct FilteredGridPartIndexSetSelector
    {
      typedef AdaptiveLeafIndexSet< FilteredGP > IndexSetType;

      static IndexSetType *create(const FilteredGP &gridPart)
      {
        return new IndexSetType( gridPart );
      }

      template < class IndexSetPtr >
      static const IndexSetType&
      indexSet ( const FilteredGP &gridPart, const std::unique_ptr< IndexSetPtr >& idxSetPtr )
      {
        assert( idxSetPtr );
        return *idxSetPtr;
      }
    };


    // FilteredGridPartIndexSetSelector
    // specialization for non-filtered index set,
    // i.e. host index set
    // -----------------------------------------

    template< class FilteredGP, class HostGP >
    struct FilteredGridPartIndexSetSelector< FilteredGP, HostGP, false >
    {
      typedef typename HostGP::IndexSetType IndexSetType;

      static IndexSetType *create(const FilteredGP &gridPart)
      {
        return nullptr;
      }

      template < class IndexSetPtr >
      static const IndexSetType&
      indexSet ( const FilteredGP &gridPart, const std::unique_ptr< IndexSetPtr >& )
      {
        return gridPart.hostGridPart().indexSet();
      }
    };



    // EntityGridTypeGetter
    // --------------------

    template< class Entity >
    struct EntityGridTypeGetter;

    template< int codim, int dim, class Grid, template< int, int, class > class Impl >
    struct EntityGridTypeGetter< Dune::Entity< codim, dim, Grid, Impl > >
    {
      typedef Grid Type;
    };

    template< class Entity >
    struct EntityGridTypeGetter< const Entity >
    {
      typedef typename EntityGridTypeGetter< Entity >::Type Type;
    };



    // FilteredGridPartTraits
    // ----------------------

    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
    struct FilteredGridPartTraits
    {
      //! \brief type of grid part
      typedef FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > GridPartType;

      struct GridPartFamily
      {
        typedef FilterImp Filter;
        typedef HostGridPartImp HostGridPart;

        static const int dimension = HostGridPart::dimension;
        static const int dimensionworld = HostGridPart::dimensionworld;

        typedef typename HostGridPart::ctype ctype;

        typedef FilteredGridPartIntersectionIterator< const GridPartFamily > IntersectionIteratorImpl;
        typedef FilteredGridPartIntersection< Filter, typename HostGridPart::IntersectionType > IntersectionImpl;

        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImpl, IntersectionImpl > IntersectionIterator;
        typedef Dune::Intersection< const GridPartFamily, IntersectionImpl > Intersection;

        template< int codim >
        struct Codim
        {
          typedef typename HostGridPart::template Codim< codim >::GeometryType Geometry;
          typedef typename HostGridPart::template Codim< codim >::LocalGeometryType LocalGeometry;

          typedef typename HostGridPart::template Codim< codim >::EntityType Entity;
          typedef typename HostGridPart::template Codim< codim >::EntitySeedType EntitySeed;
        };
      };

      //! \brief grid part imp
      typedef HostGridPartImp HostGridPartType;

      //! \brief type of grid
      typedef typename HostGridPartType::GridType GridType;

      /** \brief The type of the corresponding TwistUtility */
      typedef MetaTwistUtility< typename HostGridPartType :: TwistUtilityType >  TwistUtilityType ;

      //! \brief export filter type
      typedef FilterImp FilterType;

      //! \brief type of entity
      typedef typename FilterType::EntityType EntityType;

      //! \brief index set use in this gridpart
      typedef FilteredGridPartIndexSetSelector< GridPartType, HostGridPartType, useFilteredIndexSet > IndexSetSelectorType;

      //! \brief index set use in this gridpart
      typedef typename IndexSetSelectorType::IndexSetType IndexSetType;

      //! \brief of host grid part intersection iterator type
      typedef typename HostGridPartType::Traits::IntersectionIteratorType HostIntersectionIteratorType;

      //! \brief type of intersection iterator
      typedef typename GridPartFamily::IntersectionIterator IntersectionIteratorType;

      //! \brief type of intersection
      typedef typename GridPartFamily::Intersection IntersectionType;

      //! \brief struct providing types of the iterators on codimension cd
      template< int codim >
      struct Codim
      {
        typedef typename HostGridPartType::template Codim< codim >::GeometryType GeometryType;
        typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometryType;

        typedef typename HostGridPartType::template Codim< codim >::EntityType EntityType;
        typedef typename HostGridPartType::template Codim< codim >::EntitySeedType EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::EntityIterator< codim, typename EntityGridTypeGetter< EntityType >::Type, FilteredGridPartIterator< codim, pitype, GridPartType > > IteratorType;
        };

        typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;
      };

      typedef typename HostGridPartType::CollectiveCommunicationType CollectiveCommunicationType;

      //! \brief maximum partition type, the index set provides indices for
      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;

      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      //! \brief is true if grid on this view only has conforming intersections
      static const bool conforming = HostGridPartType::Traits::conforming;
    };



    //***************************************************************************
    //
    // FilteredGridPart
    //
    /** @addtogroup FilterGridPart
     A FilteredGridPart is a subset of a GridPart and a GridPart itself.
     The entities that belong to the FilteredGrid are defined by a
     filter class.

     Note that codim 0 entities have a method hasBoundaryIntersection().
     In general, this method will be inconsistent with the intersections
     returned by the filtered gridpart since entities are not wrapped.
    **/

    /** @ingroup FilterGridPart
     @brief
     A FilteredGridPart allows to extract a set of entities from a grid
     satisfying a given constrainted defined through a filter class.
    **/


    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
    class FilteredGridPart
    : public GridPartInterface< FilteredGridPartTraits< HostGridPartImp, FilterImp, useFilteredIndexSet > >
    {
      // type of this
      typedef FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > ThisType;

    public:
      //- Public typedefs and enums
      //! \brief traits class
      typedef FilteredGridPartTraits< HostGridPartImp, FilterImp, useFilteredIndexSet > Traits;

      //! \brief type of filter
      typedef FilterImp FilterType;

      // type of host grid part
      typedef typename Traits::HostGridPartType HostGridPartType;

      //! \brief grid type
      typedef typename Traits::GridType GridType;

      //! \brief index set type
      typedef typename Traits::IndexSetType IndexSetType;

      //! \brief intersection iterator type
      typedef typename Traits:: IntersectionIteratorType IntersectionIteratorType;

      //! \brief intersection type
      typedef typename IntersectionIteratorType::Intersection IntersectionType;

      typedef typename Traits::CollectiveCommunicationType CollectiveCommunicationType;

      //! \brief grid part typedefs, use those of traits
      template< int codim >
      struct Codim : public Traits :: template Codim< codim >
      {};

    private:
      typedef typename Traits::IndexSetSelectorType IndexSetSelectorType;

      typedef typename Codim< 0 >::EntityType EntityType;

    public:
      //- Public methods
      //! \brief constructor
      FilteredGridPart ( HostGridPartType &hostGridPart, const FilterType &filter )
      : hostGridPart_( hostGridPart ),
        filter_( filter ),
        indexSetPtr_( IndexSetSelectorType::create( *this ) )
      {
      }

      //! \brief copy constructor
      FilteredGridPart ( const FilteredGridPart &other )
      : hostGridPart_( other.hostGridPart_ ),
        filter_( other.filter_ ),
        indexSetPtr_ ( IndexSetSelectorType::create( *this ) )
      { }

      //! \brief return const reference to underlying grid
      const GridType &grid () const
      {
        return hostGridPart().grid();
      }

      //! \brief return reference to underlying grid
      GridType &grid ()
      {
        return hostGridPart().grid();
      }

      //! \brief return index set of this grid part
      //         if IndexSetType is from host grid part the original index set is returned
      const IndexSetType &indexSet() const
      {
        return IndexSetSelectorType::indexSet( *this, indexSetPtr_ );
      }

      //! \brief Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim >::IteratorType begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      //! \brief Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType begin () const
      {
        typedef typename Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
        typedef FilteredGridPartIterator< codim, pitype, ThisType > IteratorImpl;
        return IteratorType( IteratorImpl( *this, hostGridPart().template begin< codim, pitype >() ) );
      }

      //! \brief Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim >::IteratorType end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! \brief End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType end () const
      {
        typedef typename Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
        typedef FilteredGridPartIterator< codim, pitype, ThisType > IteratorImpl;
        return IteratorType( IteratorImpl( *this, hostGridPart().template end< codim, pitype >() ) );
      }

      //! \brief Returns maxlevel of the grid
      int level () const
      {
        return hostGridPart().level();
      }

      //! \brief ibegin of corresponding intersection iterator for given entity
      IntersectionIteratorType ibegin ( const EntityType &entity ) const
      {
        typedef typename IntersectionIteratorType::Implementation IntersectionIteratorImpl;
        return IntersectionIteratorType( IntersectionIteratorImpl( filter(), hostGridPart().ibegin( entity ) ) );
      }

      //! \brief iend of corresponding intersection iterator for given entity
      IntersectionIteratorType iend ( const EntityType &entity ) const
      {
        typedef typename IntersectionIteratorType::Implementation IntersectionIteratorImpl;
        return IntersectionIteratorType( IntersectionIteratorImpl( filter(), hostGridPart().iend( entity ) ) );
      }

      //! \brief boundary id
      int boundaryId ( const IntersectionType &intersection ) const
      {
        return hostGridPart().boundaryId( intersection.impl().hostIntersection() );
      }

      const CollectiveCommunicationType &comm () const { return hostGridPart_.comm(); }

      //! \brief corresponding communication method for this grid part
      template < class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &dataHandle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandleImp, DataType >  HostHandleType;
        FilteredGridPartDataHandle< HostHandleType, ThisType > handleWrapper( dataHandle, *this );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      /** \copydoc GridPartInterface::entity(const EntitySeed &seed) const */
      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return hostGridPart().entity( seed );
      }

      //! \brief return reference to filter
      const FilterType &filter () const
      {
        return filter_;
      }

      //! \brief return reference to filter
      FilterType &filter ()
      {
        return filter_;
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return filter().contains( entity );
      }

      HostGridPartType &hostGridPart ()
      {
        return hostGridPart_;
      }

      const HostGridPartType &hostGridPart () const
      {
        return hostGridPart_;
      }

      /** \copydoc GridPartInterface::convert(const Entity &entity) const */
      template <class Entity>
      const Entity& convert ( const Entity &entity ) const
      {
        return hostGridPart().convert( entity );
      }

    private:
      HostGridPartType &hostGridPart_;
      FilterType filter_;
      std::unique_ptr< IndexSetType > indexSetPtr_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_HH
