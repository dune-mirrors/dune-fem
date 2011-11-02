#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_DEADITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_DEADITERATOR_HH

#include <dune/grid/common/entityiterator.hh>

#include <dune/fem/gridpart/idgridpart/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // DeadIterator
    // ------------
    
    template< int codim, class GridFamily >
    class DeadIterator
    : public IdEntityPointer< IdEntityPointerTraits< codim, GridFamily > >
    {
      typedef IdEntityPointer< IdEntityPointerTraits< codim, GridFamily > > Base;

    protected:
      typedef typename Base::HostIteratorType HostIteratorType;

      using Base::hostIterator_;
      using Base::releaseEntity;

    public:
      explicit DeadIterator ( const HostIteratorType &hostIterator )
      : Base( hostIterator )
      {}
      
      void increment ()
      {
        DUNE_THROW( InvalidStateException, "Trying to increment a dead iterator." );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_DEADITERATOR_HH
