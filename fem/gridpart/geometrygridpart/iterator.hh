#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH

#include <cassert>

#include <type_traits>
#include <utility>

#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartIterator
    // -----------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    class GeometryGridPartIterator
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

    public:
      typedef typename Traits::GridFunctionType GridFunctionType;
      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;

      static const int codimension = HostIteratorType::codimension;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

      GeometryGridPartIterator () = default;

      GeometryGridPartIterator ( const GridFunctionType &gridFunction, HostIteratorType hostIterator )
      : gridFunction_( &gridFunction ),
        hostIterator_( std::move( hostIterator ) )
      {}

      void increment ()
      {
        ++hostIterator_;
      }

      Entity dereference () const
      {
        return typename Entity::Implementation( gridFunction(), *hostIterator_ );
      }

      bool equals ( const GeometryGridPartIterator &rhs ) const
      {
        return hostIterator_ == rhs.hostIterator_;
      }

      int level () const
      {
        return hostIterator_.level();
      }

      operator Dune::DefaultEntityPointer< Entity > () const
      {
        return Dune::DefaultEntityPointer< Entity >( dereference() );
      }

      bool equals ( const Dune::DefaultEntityPointer< Entity > &rhs ) const
      {
        return dereference() == rhs.dereference();
      }

    private:
      const GridFunctionType &gridFunction () const
      {
        assert( gridFunction_ );
        return *gridFunction_;
      }

      const GridFunctionType *gridFunction_ = nullptr;
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH
