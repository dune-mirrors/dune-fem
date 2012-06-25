#include <config.h>

#include "basefunctions.hh"

namespace Dune
{
  namespace Fem
  {

  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 1 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 1 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 1 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 1 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 1 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 1 >;

  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 2 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 2 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 2 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 2 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 2 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 2 >;

  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 3 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 3 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 3 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 3 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 3 >;
  template class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 3 >;

  } // namespace Fem
  
} // namespace Dune
