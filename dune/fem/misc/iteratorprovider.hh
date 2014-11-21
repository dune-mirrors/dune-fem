#ifndef DUNE_FEM_MISC_ITERATORPROVIDER_HH
#define DUNE_FEM_MISC_ITERATORPROVIDER_HH

#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  namespace Fem
  {

    // IteratorProvider
    // ----------------

    template< class DiscreteFunctionSpace >
    struct IteratorProvider
    {
      typedef typename DiscreteFunctionSpace::IteratorType IteratorType;

      explicit IteratorProvider ( const DiscreteFunctionSpace &space )
      : space_( space )
      {}

      IteratorType begin () const { return space_.begin(); }
      IteratorType end () const { return space_.end(); }

    private:
      const DiscreteFunctionSpace &space_;
    };


    // PartitionIteratorProvider
    // -------------------------

    template< class DiscreteFunctionSpace, PartitionIteratorType pitype >
    struct PartitionIteratorProvider
    {
      typedef typename DiscreteFunctionSpace::GridPartType GridPartType;

      static const int codimension = DiscreteFunctionSpace::Traits::codimension;
      typedef typename GridPartType::template Codim< codimension >::template Partition< pitype >::IteratorType IteratorType;

      explicit PartitionIteratorProvider ( const DiscreteFunctionSpace &space )
      : gridPart_( space.gridPart() )
      {}

      IteratorType begin () const { return gridPart_.template begin< codimension, pitype >(); }
      IteratorType end () const { return gridPart_.template end< codimension, pitype >(); }

    private:
      const GridPartType &gridPart_;
    };
  } // namespace Fem

} // namespace Dune

#endif  // #ifndef DUNE_FEM_MISC_ITERATORPROVIDER_HH
