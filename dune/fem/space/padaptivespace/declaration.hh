#ifndef DUNE_FEM_SPACE_PADAPTIVE_DECLARATION_HH
#define DUNE_FEM_SPACE_PADAPTIVE_DECLARATION_HH

namespace Dune
{

  namespace Fem 
  {

    // PAdaptiveDGSpace
    // ----------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storagee >
    class PAdaptiveDGSpace;


    // PAdaptiveLagrangeSpace
    // ----------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class PAdaptiveLagrangeSpace;

  } // namespace Fem 

}  // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_DECLARATION_HH
