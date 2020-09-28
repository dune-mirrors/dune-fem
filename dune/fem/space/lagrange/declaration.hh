#ifndef DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH
#define DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH

#include <dune/fem/space/shapefunctionset/selectcaching.hh>

namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscreteFunctionSpace
    // -----------------------------

    template< class FunctionSpace, class GridPart, int order, class Storage >
    class LagrangeDiscreteFunctionSpace;

    // 6 is the maximal polynomial order that can be selected with DynamicLagrangeDiscreteFunctionSpace
    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    using DynamicLagrangeDiscreteFunctionSpace = LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, -6, Storage >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH
