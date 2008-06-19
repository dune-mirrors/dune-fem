#include <config.h>

#include "orthonormalbase_1d.cc"
#include "dgbasefunctions.hh"

namespace Dune
{

  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 0 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 1 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 2 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 3 >;

}
