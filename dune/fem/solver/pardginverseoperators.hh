#ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH

#include <cassert>

#include <limits>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>

namespace Dune
{

  namespace Fem
  {

    // ParDGGeneralizedMinResInverseOperator
    // -------------------------------------
    template< class DiscreteFunction >
    using  ParDGGeneralizedMinResInverseOperator = GmresInverseOperator< DiscreteFunction >;

    // ParDGBiCGStabInverseOperator
    // ----------------------------
    template< class DiscreteFunction >
    using ParDGBiCGStabInverseOperator = BicgstabInverseOperator< DiscreteFunction >;
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
