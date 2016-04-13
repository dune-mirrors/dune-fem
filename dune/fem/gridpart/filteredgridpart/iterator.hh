#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ITERATOR_HH

//- system includes
#include <cassert>

//- dune-grid includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartIterator
    // ------------------------

    template< int codim, PartitionIteratorType pitype, class GridPartImp >
    class FilteredGridPartIterator
    {
      typedef FilteredGridPartIterator< codim, pitype, GridPartImp > ThisType;

      typedef GridPartImp GridPartType;
      typedef typename GridPartType::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;

    public:
      // type of entity
      typedef typename HostIteratorType::Entity Entity;

      static const int codimension = codim;

      //! \brief constructor
      FilteredGridPartIterator ( const GridPartType &gridPart, const HostIteratorType &hostIterator )
      : gridPart_( gridPart ),
        hostIterator_( hostIterator ),
        hostEnd_( gridPart.hostGridPart().template end< codim, pitype >() )
      {
        if( done() )
          return;

        if( !gridPart.contains( *hostIterator_ ) )
          increment();
      }

      //! \brief constructor
      FilteredGridPartIterator ( const ThisType &other )
      : gridPart_( other.gridPart_ ),
        hostIterator_( other.hostIterator_ ),
        hostEnd_( other.hostEnd_ )
      {}

      //! \brief assignment operator
      ThisType &operator= ( const ThisType &other )
      {
        assert( &gridPart_ == &other.gridPart_ );
        hostIterator_ = other.hostIterator_;
        hostEnd_ = other.hostEnd_;
        return *this;
      }

      //! \brief increment
      void increment ()
      {
        assert( !done() );
        do { ++hostIterator_; } while ( !done() && !contains() );
      }

      //! \brief return level
      int level () const { return hostIterator_.level(); }

      //! return reference to entity object
      Entity dereference () const { return *hostIterator_; }

      //! \brief check for equality
      bool equals ( const ThisType &other ) const
      {
        return hostIterator_ == other.hostIterator_;
      }

    private:
      bool done () const
      {
        return (hostIterator_ == hostEnd_);
      }

      bool contains () const
      {
        assert( !done() );
        return gridPart().contains( *hostIterator_ );
      }

      // reference to grid part
      const GridPartType &gridPart () const { return gridPart_; }

      const GridPartType &gridPart_;
      HostIteratorType hostIterator_;
      HostIteratorType hostEnd_;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ITERATOR_HH
