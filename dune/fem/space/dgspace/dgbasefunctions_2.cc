#include <config.h>

#include "orthonormalbase_2d.cc"
#include "dgbasefunctions.hh"

namespace Dune
{

  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 0 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 1 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 2 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 3 >;
  
}
