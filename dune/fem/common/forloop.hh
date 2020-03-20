#ifndef DUNE_FEM_FORLOOP_HH
#define DUNE_FEM_FORLOOP_HH

#include <utility>

#include <dune/common/deprecated.hh>
#include <dune/common/hybridutilities.hh>

namespace Dune
{

  namespace Fem {

  template< template< int > class Operation, int first, int last >
  struct ForLoop
  {
    static_assert( (first <= last), "Fem::Fem::ForLoop: first > last" );

    template<typename... Args>
    static void apply(Args&&... args)
    {
      Dune::Hybrid::forEach(std::make_index_sequence<last+1-first>{},
        [&](auto i){Operation<i+first>::apply(args...);});
    }
  };

  } // end namespace Fem
} // end namespace Dune
#endif // #ifndef DUNE_FEM_FORLOOP_HH
