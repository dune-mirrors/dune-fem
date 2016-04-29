#ifndef DUNE_FEM_IO_PARAMETER_EXCEPTIONS_HH
#define DUNE_FEM_IO_PARAMETER_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

namespace Dune
{

  namespace Fem
  {

    // ParameterNotFound
    // -----------------

    class ParameterNotFound
      : public Exception
    {};



    // ParameterInvalid
    // ----------------

    class ParameterInvalid
      : public Exception
    {};

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_EXCEPTIONS_HH
