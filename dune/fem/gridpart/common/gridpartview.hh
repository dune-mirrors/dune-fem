#ifndef DUNE_FEM_GRIDVIEW_HH
#define DUNE_FEM_GRIDVIEW_HH

#include <dune/grid/common/gridview.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart >
    class GridPartViewImpl;


    template< class GridPart >
    struct GridPartViewTraits
    {
      typedef GridPartViewImpl< GridPart > GridViewImp;

      typedef typename GridPart::GridType Grid;
      typedef typename GridPart::IndexSetType IndexSet;
      typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

      typedef typename IntersectionIterator::Intersection Intersection;

      typedef typename GridPart::CollectiveCommunicationType CollectiveCommunication;

      template< int codim >
      struct Codim
      : public Grid::Traits::template Codim< codim >
      {
        typedef typename GridPart::template Codim< codim >::IteratorType Iterator;

        typedef typename GridPart::template Codim< codim >::EntityType Entity;
        typedef typename GridPart::template Codim< codim >::EntityPointerType EntityPointer;

        typedef typename GridPart::template Codim< codim >::GeometryType Geometry;
        typedef typename GridPart::template Codim< codim >::LocalGeometryType LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridPart::template Codim< codim >::template Partition< pitype >::IteratorType Iterator;
        };
      };

      static const bool conforming = GridPart::Traits::conforming;
    };


    template< class GridPart >
    class GridPartViewImpl
    {
      typedef GridPartViewImpl< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      typedef GridPartViewTraits< GridPartType > Traits;

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

      explicit GridPartViewImpl ( const GridPartType &gridPart )
      : gridPart_( gridPart )
      {}

      GridPartViewImpl ( const ThisType &other )
      : gridPart_( other.gridPart_ )
      {}

    private:
      ThisType &operator= ( const ThisType & );

    public:
      const Grid &grid () const
      {
        return gridPart_.grid();
      }

      const IndexSet &indexSet () const
      {
        return gridPart_.indexSet();
      }

      int size ( int codim ) const
      {
        return indexSet().size( codim );
      }

      int size ( const GeometryType &type ) const
      {
        return indexSet().size( type );
      }

      template< int codim >
      typename Codim< codim >::Iterator begin () const
      {
        return gridPart_.template begin< codim >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator begin () const
      {
        return gridPart_.template begin< codim, pitype >();
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        return gridPart_.template end< codim >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator end () const
      {
        return gridPart_.template end< codim, pitype >();
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        return gridPart_.ibegin( entity );
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        return gridPart_.iend( entity );
      }

      const CollectiveCommunication &comm () const
      {
        return gridPart_.comm();
      }

      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype,
                         CommunicationDirection dir ) const
      {
        gridPart_.communicate( data, iftype, dir );
      }

    private:
      const GridPartType &gridPart_;
    };



    template< class GridPart >
    class GridPartView
    : public GridView< GridPartViewTraits< GridPart > >
    {
      typedef GridPartView< GridPart > ThisType;
      typedef GridView< GridPartViewTraits< GridPart > > BaseType;

      typedef typename BaseType::GridViewImp GridViewImp;

    public:
      explicit GridPartView ( const GridPart &gridPart )
      : BaseType( GridViewImp( gridPart ) )
      {}

      GridPartView ( const ThisType &other )
      : BaseType( other )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDVIEW_HH
