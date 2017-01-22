#ifndef DUNE_FEM_GRIDPART_COMMON_GRIDPART2GRIDVIEW_HH
#define DUNE_FEM_GRIDPART_COMMON_GRIDPART2GRIDVIEW_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart >
    class GridPart2GridViewImpl;


    template< class GridPart >
    struct GridPart2GridViewTraits
    {
      typedef GridPart2GridViewImpl< GridPart > GridViewImp;

      typedef typename GridPart::GridType Grid;
      typedef typename GridPart::IndexSetType IndexSet;
      typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

      typedef typename IntersectionIterator::Intersection Intersection;

      typedef typename GridPart::CollectiveCommunicationType CollectiveCommunication;

      template< int codim >
      struct Codim
      : public Grid::Traits::template Codim< codim >
      {
        typedef typename GridPart::template Codim< codim >::EntityType Entity;

        typedef typename GridPart::template Codim< codim >::GeometryType Geometry;
        typedef typename GridPart::template Codim< codim >::LocalGeometryType LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridPart::template Codim< codim >::template Partition< pitype >::IteratorType Iterator;
        };

        typedef typename Partition< All_Partition >::Iterator Iterator;
      };

      static const bool conforming = GridPart::Traits::conforming;
    };


    template< class GridPart >
    class GridPart2GridViewImpl
    {
      typedef GridPart2GridViewImpl< GridPart > ThisType;

    public:
      typedef typename GridPart::ctype ctype;

      typedef GridPart GridPartType;

      typedef GridPart2GridViewTraits< GridPartType > Traits;

      /** \brief type of the grid */
      typedef typename Traits::Grid Grid;

      /** \brief type of the index set */
      typedef typename Traits::IndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Traits::Intersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Traits::IntersectionIterator IntersectionIterator;

      /** \brief type of the collective communication */
      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      /** \brief Codim Structure */
      template< int codim >
      struct Codim
      : public Traits::template Codim< codim >
      {};

      enum { conforming = Traits::conforming };

      enum { dimension = GridPartType::dimension };
      enum { dimensionworld = GridPartType::dimensionworld };

      explicit GridPart2GridViewImpl ( const GridPartType &gridPart )
        : gridPart_( &gridPart )
      {}

      const Grid &grid () const
      {
        return gridPart().grid();
      }

      const IndexSet &indexSet () const
      {
        return gridPart().indexSet();
      }

      int size ( int codim ) const
      {
        return indexSet().size( codim );
      }

      int size ( const GeometryType &type ) const
      {
        return indexSet().size( type );
      }

      template<class EntityType>
      bool contains (const EntityType& e) const
      {
        return indexSet().contains(e);
      }

      template< int codim >
      typename Codim< codim >::Iterator begin () const
      {
        return begin< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator begin () const
      {
        return gridPart().template begin< codim, pitype >();
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        return end< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator end () const
      {
        return gridPart().template end< codim, pitype >();
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        return gridPart().ibegin( entity );
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        return gridPart().iend( entity );
      }

      const CollectiveCommunication &comm () const
      {
        return gridPart().comm();
      }

      int overlapSize ( int codim ) const
      {
        DUNE_THROW( NotImplemented, "Method ghostSize() not implemented yet" );
      }

      int ghostSize( int codim ) const
      {
        DUNE_THROW( NotImplemented, "Method ghostSize() not implemented yet" );
      }

      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype,
                         CommunicationDirection dir ) const
      {
        gridPart().communicate( data, iftype, dir );
      }

      const GridPartType &gridPart () const { assert( gridPart_ ); return *gridPart_; }

    private:
      const GridPartType *gridPart_;
    };



    template< class GridPart >
    class GridPart2GridView
    : public GridView< GridPart2GridViewTraits< GridPart > >
    {
      typedef GridPart2GridView< GridPart > ThisType;
      typedef GridView< GridPart2GridViewTraits< GridPart > > BaseType;

      typedef typename BaseType::GridViewImp GridViewImp;

    public:
      explicit GridPart2GridView ( const GridPart &gridPart )
      : BaseType( GridViewImp( gridPart ) )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_GRIDPART2GRIDVIEW_HH
