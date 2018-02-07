#ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
#define DUNE_FEM_DEFAULTQUADRATURES_HH

//#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/quadrature/idprovider.hh>

// don't use quadratures from dune-grid
//#define USE_DUNE_QUADRATURES

// include quadrature points
#ifdef USE_DUNE_QUADRATURES
# warning "Using Dune::Geomtry Quadratures by define of USE_DUNE_QUADRATURES!"
# include "dunequadratures.hh"
#else
# include "femquadratures.hh"
#endif

namespace Dune
{

  namespace Fem
  {

#ifdef USE_DUNE_QUADRATURES
    template<class Field, int dim>
    using DefaultQuadratureTraits = DuneQuadratureTraits<Field, dim>;
#else
    template<class Field, int dim>
    using DefaultQuadratureTraits = FemQuadratureTraits<Field, dim>;
#endif

  } // namespace Fem

} // namespace Dune

#undef USE_DUNE_QUADRATURES
#endif // #ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
