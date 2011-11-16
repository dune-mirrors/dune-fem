#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH

#include <dune/grid/common/entityiterator.hh>
#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/referenceelements.hh>
#else
#include <dune/grid/common/genericreferenceelements.hh>
#endif

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
      typedef IdEntityPointer< IdIteratorTraits< codim, pitype, GridFamily > > Base;

    protected:
      typedef typename Base::HostIteratorType HostIteratorType;

      using Base::hostIterator_;
      using Base::releaseEntity;

    public:
      IdIterator ( const HostIteratorType &hostIterator )
      : Base( hostIterator )
      {}
      
      void increment ()
      {
        ++hostIterator_;
        releaseEntity();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ITERATOR_HH
