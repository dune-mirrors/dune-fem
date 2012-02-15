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

    template< int codim, PartitionIteratorType pitype, class GridPartImp, class HostIteratorImp >
    class FilteredGridPartIterator
    {
      // type of this
      typedef FilteredGridPartIterator< codim, pitype, GridPartImp, HostIteratorImp > ThisType;

      // grid part type
      typedef GridPartImp GridPartType;

      // host iterator type
      typedef HostIteratorImp HostIteratorType;

      // entity pointer type
      typedef typename GridPartType::GridType::template Codim< codim >::EntityPointer EntityPointerType;

    public:
      // type of entity
      typedef typename HostIteratorType::Entity Entity;

      //! \brief constructor
      FilteredGridPartIterator ( const GridPartType & gridPart, const HostIteratorType & hostIterator )
      : gridPart_( gridPart ),
        hostIterator_( hostIterator ),
        hostEnd_( gridPart.hostGridPart().template end< codim, pitype >() )
      {
        if( done() ) 
          return;

        if( !gridPart.contains( *hostIterator_ ) )
          ++(*this);
      }

      //! \brief constructor
      FilteredGridPartIterator ( const ThisType & other )
      : gridPart_( other.gridPart_ ),
        hostIterator_( other.hostIterator_ ),
        hostEnd_( other.hostEnd_ )
      { }

      //! \brief assignment operator
      ThisType & operator= ( const ThisType & other )
      {
        assert( &gridPart_ == &other.gridPart_ );
        hostIterator_ = other.hostIterator_;
        hostEnd_ = other.hostEnd_;
        return *this;
      }

      //! \brief increment
      ThisType & operator++ ()
      {
        assert( !done() );
        do
        {
          ++hostIterator_;
        } while ( !done() && !contains() );
        return *this;
      }

      //! \brief return level
      int level () const
      {
        return hostIterator_.level();
      }

      const Entity & operator* () const
      {
        return *hostIterator_;
      }

      const Entity * operator-> () const
      {
        return &(*hostIterator_);
      }

      //! \brief cast to entity pointer
      operator EntityPointerType & ()
      {
        return hostIterator_;
      }

      //! \brief cast to const entity pointer
      operator const EntityPointerType & () const
      {
        return hostIterator_;
      }

      //! \brief check for equality
      bool operator== ( const ThisType & other ) const
      {
        return hostIterator_.operator==( other.hostIterator_ );
      }

      //! \brief check for inequality
      bool operator != ( const ThisType & other ) const
      {
        return !(*(this)==other);
      }

    private:
      // return true for end iterator
      bool done () const
      {
        return (hostIterator_ == hostEnd_ );
      }

      bool contains () const
      {
        assert( !done() );
        return gridPart().contains( *hostIterator_ );
      }

      // reference to grid part
      const GridPartType & gridPart () const
      {
        return gridPart_;
      }

      const GridPartType & gridPart_;
      HostIteratorType hostIterator_;
      HostIteratorType hostEnd_;

    }; // end class FilteredGridPartIterator

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_ITERATOR_HH
