#include <config.h>

#include "orthonormalbase_3d.cc"
#include "dgbasefunctions.hh"

namespace Dune
{

  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 0 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 1 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 2 >;
  template class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 3 >;

}
