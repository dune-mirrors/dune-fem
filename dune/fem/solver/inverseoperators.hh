#ifndef DUNE_FEM_INVERSEOPERATORS_HH
#define DUNE_FEM_INVERSEOPERATORS_HH

// deprecated header...

#include <dune/fem/solver/cginverseoperator.hh>

namespace Dune {

#ifdef DUNE_FEM_COMPATIBILITY
// to be removed after release of version 1.3 
#define CGInverseOp  CGInverseOperator
#endif

} // namespace Dune

#endif // #ifndef DUNE_FEM_INVERSEOPERATORS_HH
