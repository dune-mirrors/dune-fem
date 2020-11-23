#ifndef DUNE_FEM_FORLOOP_HH
#define DUNE_FEM_FORLOOP_HH

#include <utility>

#include <dune/common/visibility.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/hybridutilities.hh>

namespace Dune
{

  namespace Fem {

  template< template< int > class Operation, int first, int last >
  struct ForLoop
  {
    static_assert( (first <= last), "Fem::Fem::ForLoop: first > last" );

    static const std::size_t N = last + 1 - first;

    template<typename... Args>
    static DUNE_PRIVATE void apply(Args&&... args)
    {
      Dune::Hybrid::forEach(std::make_index_sequence<N>{},
        [&](auto i){Operation<i+first>::apply(args...);});
    }
  };

  } // end namespace Fem
} // end namespace Dune
#endif // #ifndef DUNE_FEM_FORLOOP_HH
