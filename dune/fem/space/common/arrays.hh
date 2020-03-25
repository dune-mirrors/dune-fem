#ifndef DUNE_FEM_ARRAYS_HH
#define DUNE_FEM_ARRAYS_HH

#warning("WARNING: MutableArray defined in dune/fem/space/common/arrays.hh is deprecated. Please use directly DynamicArray defined in dune/fem/storage/dynamicarray.hh.")

#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{
  namespace Fem
  {
    template<typename... Args>
    using MutableArray = DynamicArray<Args...>;
  }
}

#endif // DUNE_FEM_ARRAYS_HH
