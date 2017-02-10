#ifndef DUNE_FEM_SPACE_FOURIER_DECLARATION_HH
#define DUNE_FEM_SPACE_FOURIER_DECLARATION_HH

#include <type_traits>

namespace Dune
{

  namespace Fem
  {

    // FourierDiscreteFunctionSpace
    // ----------------------------

    template< class FunctionSpace, class GridPart, int order >
    class FourierDiscreteFunctionSpace;



    // IsFourierDiscreteFunctionSpace
    // ------------------------------

    template< class DiscreteFunctionSpace >
    struct IsFourierDiscreteFunctionSpace
      : std::integral_constant< bool, false >
    {};

    template< class FunctionSpace, class GridPart, int order >
    struct IsFourierDiscreteFunctionSpace< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, order > >
      : std::integral_constant< bool, true >
    {};

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_DECLARATION_HH
