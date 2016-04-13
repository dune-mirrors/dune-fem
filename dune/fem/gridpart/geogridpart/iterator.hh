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

    // GeoIterator
    // -----------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    class GeoIterator
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

    public:
      typedef typename Traits::CoordFunctionType CoordFunctionType;
      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;

      static const int codimension = HostIteratorType::codimension;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

      GeoIterator () = default;

      GeoIterator ( const CoordFunctionType &coordFunction, HostIteratorType hostIterator )
      : coordFunction_( &coordFunction ),
        hostIterator_( std::move( hostIterator ) )
      {}

      void increment ()
      {
        ++hostIterator_;
      }

      Entity dereference () const
      {
        return typename Entity::Implementation( coordFunction(), *hostIterator_ );
      }

      bool equals ( const GeoIterator &rhs ) const
      {
        return hostIterator_ == rhs.hostIterator_;
      }

      int level () const
      {
        return hostIterator_.level();
      }

    private:
      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const CoordFunctionType *coordFunction_ = nullptr;
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH
