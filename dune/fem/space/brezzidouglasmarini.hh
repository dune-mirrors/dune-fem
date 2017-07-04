#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH

#ifdef USE_OLD_BDM_SPACE

#include "brezzidouglasmarini/space.hh"

#else

#include <dune/fem/space/localfiniteelement/space.hh>
#include <dune/fem/space/brezzidouglasmarini/localfemap.hh>

namespace Dune
{
  namespace Fem
  {

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    using BDMDiscreteFunctionSpace
    = LocalFiniteElementDiscreteFunctionSpace< BDMLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >,
      FunctionSpace, Storage >;
  }
}
#endif
#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
