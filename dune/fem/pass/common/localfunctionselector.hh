#ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONSELECTOR_HH
#define DUNE_FEM_PASS_COMMON_LOCALFUNCTIONSELECTOR_HH

#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/function/localfunction/const.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionSelector
    // ---------------------

    /*
     * \brief Please doc me.
     */
    template< class DiscreteFunctionType >
    struct LocalFunctionSelector
    {
      typedef ConstLocalFunction< DiscreteFunctionType > Type;
    };

    template< class Evaluate >
    struct LocalFunctionSelector< LocalFunctionAdapter< Evaluate > >
    {
      typedef typename LocalFunctionAdapter< Evaluate >::LocalFunctionType Type;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONSELECTOR_HH
