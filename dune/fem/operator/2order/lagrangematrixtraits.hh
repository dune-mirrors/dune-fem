#ifndef DUNE_FEM_LAGRANGEMATRIXTRAITS_HH
#define DUNE_FEM_LAGRANGEMATRIXTRAITS_HH

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  namespace Fem
  {

    template< class RowSpace, class ColumnSpace = RowSpace,
              bool addNonConformingNeighbors = false >
    struct LagrangeMatrixTraits
    {
      typedef RowSpace RowSpaceType;
      typedef ColumnSpace ColumnSpaceType;
      typedef Dune::LagrangeMatrixSetup< addNonConformingNeighbors > StencilType;
      typedef Dune::ParallelScalarProduct< ColumnSpaceType > ParallelScalarProductType;

#if HAVE_DUNE_ISTL
      template< class M >
      struct Adapter
      {
        typedef Dune::LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
      };
#endif // #HAVE_DUNE_ISTL
    };

  } // namespace Dune

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGEMATRIXTRAITS_HH
