#include <config.h>

#include "orthonormalbase_mod.cc"
#include "dgbasefunctions.hh"

namespace Dune
{

  namespace Library
  {

    template< unsigned int dim, int pOrder >
    void Instantiate_DGBaseFunctionFactory ()
    {
      DiscontinuousGalerkinBaseFunctionFactory
        < FunctionSpace< double, double, dim, 1 >, pOrder >
        factory( GeometryType( GeometryType :: simplex, dim ) );
      factory.numBaseFunctions();
    }
    
    void Instantiate_DGSpace_Templates ()
    {
      Instantiate_DGBaseFunctionFactory< 2, 0 >();
      Instantiate_DGBaseFunctionFactory< 2, 1 >();
      Instantiate_DGBaseFunctionFactory< 2, 2 >();
      Instantiate_DGBaseFunctionFactory< 2, 3 >();
    }
  }
  
}
