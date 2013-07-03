#ifndef DUNE_FEM_DGMATRIXTRAITS_HH
#define DUNE_FEM_DGMATRIXTRAITS_HH

#warning DEPRECATED file not needed anymore - use new XXXLinearOperator instead of XXXMatrixOperator

#include <dune/fem/operator/2order/dgmatrixsetup.hh>

namespace Dune
{

  namespace Fem
  {
    template< class RowSpace, class ColumnSpace = RowSpace >
    struct DGMatrixTraits
    {
      typedef RowSpace RowSpaceType;
      typedef ColumnSpace ColumnSpaceType;
      typedef Dune::ElementAndNeighbors StencilType;
      typedef ParallelScalarProduct< ColumnSpaceType > ParallelScalarProductType;

#if HAVE_DUNE_ISTL
      template< class M >
      struct Adapter
      {
        typedef Dune::DGParallelMatrixAdapter< M > MatrixAdapterType;
      };
#endif // #HAVE_DUNE_ISTL
    };

  } // namespace Dune

} // namespace Dune

#endif // #ifndef DUNE_FEM_DGMATRIXTRAITS_HH
