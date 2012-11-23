#ifndef DUNE_FEM_SPACE_GENERICDISCRETE_DECLARATION_HH
#define DUNE_FEM_SPACE_GENERICDISCRETE_DECLARATION_HH

namespace Dune
{

  namespace Fem
  {

    // GenericDiscreteFunctionSpace
    // -----------------------------------

    template< class GridPart, class LocalFiniteElement, int polOrder, template< class > class Storage >
    struct GenericDiscreteFunctionSpace;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_GENERICDISCRETE_DECLARATION_HH
