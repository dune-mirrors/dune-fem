#ifndef DUNE_FEMPY_GRID_ADAPTIVEDOFVECTOR_HH
#define DUNE_FEMPY_GRID_ADAPTIVEDOFVECTOR_HH

#include <cstddef>

#include <memory>
#include <utility>

#include <dune/fempy/grid/virtualizedrestrictprolong.hh>

namespace Dune
{

  namespace FemPy
  {

    // AdaptiveDofVector
    // -----------------

    template< class Grid, class D = double >
    struct AdaptiveDofVector
    {
      typedef D DofType;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      typedef VirtualizedRestrictProlong< Grid > RestrictProlong;

      virtual ~AdaptiveDofVector () = default;

      virtual void enableDofCompression () = 0;

      virtual std::size_t numLocalDofs ( const Element &element ) const = 0;

      virtual DofType *getLocalDofs ( const Element &element, DofType *localDofs ) const = 0;
      virtual const DofType *setLocalDofs ( const Element &element, const DofType *localDofs ) = 0;

      virtual RestrictProlong restrictProlong () = 0;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_ADAPTIVEDOFVECTOR_HH
