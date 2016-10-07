#ifndef DUNE_FEM_THREADFILTER_HH
#define DUNE_FEM_THREADFILTER_HH

#warning("WARNING: please use DomainFilter instead of ThreadFilter.")

#include <dune/fem/gridpart/filter/domainfilter.hh>

namespace Dune
{
  namespace Fem
  {
    template< typename ... Args >
    using ThreadFilter = DomainFilter< Args ... >;
  }
}

#endif // #ifndef DUNE_FEM_THREADFILTER_HH
