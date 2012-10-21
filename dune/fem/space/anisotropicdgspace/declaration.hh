#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DECLARATION_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DECLARATION_HH

/**
  @file
  @author Christoph Gersbacher
  @brief Provides a forward declaration of anisotropic DG space 
*/


namespace Dune
{

  namespace Fem 
  {

    // AnisotropicDGSpace
    // ------------------

    template< class FunctionSpace, class GridPart, int maxOrder, template< class > class Storage >
    class AnisotropicDGSpace;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DECLARATION_HH
