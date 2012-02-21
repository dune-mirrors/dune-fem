#ifndef DUNE_FEM_SPACE_DOFMAPPER_EXCEPTIONS_HH
#define DUNE_FEM_SPACE_DOFMAPPER_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

namespace Dune
{

  namespace Fem
  {

    class DofMapperError
    : public Dune::Exception
    {};

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DOFMAPPER_EXCEPTIONS_HH
