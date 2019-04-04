#ifndef DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH
#define DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH

namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscreteFunctionSpace
    // -----------------------------

    template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
    class LagrangeDiscreteFunctionSpace;

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    using DynamicLagrangeDiscreteFunctionSpace = LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, -1, Storage >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_DECLARATION_HH
