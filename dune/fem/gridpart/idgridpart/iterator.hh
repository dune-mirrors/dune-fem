#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH

#include <dune/grid/common/entityiterator.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/fem/gridpart/idgridpart/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIteratorTraits
    // ----------------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    struct IdIteratorTraits
    : public IdEntityPointerTraits< codim, GridFamily >
    {
      typedef typename remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;
    };



    // IdIterator
    // ----------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    class IdIterator
    : public IdEntityPointer< IdIteratorTraits< codim, pitype, GridFamily > >
    {
      typedef IdEntityPointer< IdIteratorTraits< codim, pitype, GridFamily > > BaseType;

    protected:
      typedef typename BaseType::HostIteratorType HostIteratorType;
      typedef typename BaseType::ExtraData        ExtraData;

      using BaseType::hostIterator_;

    public:
      IdIterator ( ExtraData data, HostIteratorType hostIterator )
      : BaseType( data, hostIterator )
      {}

      void increment ()
      {
        ++hostIterator_;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
