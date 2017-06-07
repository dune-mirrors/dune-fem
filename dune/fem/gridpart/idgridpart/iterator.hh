#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/version.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIterator
    // ----------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    class IdIterator
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

    public:
      typedef typename Traits::ExtraData ExtraData;
      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;

#if !DUNE_VERSION_NEWER(DUNE_GRID, 2, 6 )
      static const int codimension = HostIteratorType::codimension;
#endif

      typedef typename Traits::template Codim< codim >::Entity Entity;

      IdIterator () = default;

      IdIterator ( ExtraData data, HostIteratorType hostIterator )
      : data_( std::move( data ) ),
        hostIterator_( std::move( hostIterator ) )
      {}

      void increment ()
      {
        ++hostIterator_;
      }

      Entity dereference () const
      {
        return typename Entity::Implementation( data_, *hostIterator_ );
      }

      bool equals ( const IdIterator &rhs ) const
      {
        return hostIterator_ == rhs.hostIterator_;
      }

      int level () const
      {
        return hostIterator_.level();
      }

    private:
      ExtraData data_;
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
