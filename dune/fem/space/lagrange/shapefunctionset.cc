#include <config.h>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>

namespace Dune
{

  namespace Fem
  {

    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 1 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 1 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 1 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 1 >;

    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 2 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 2 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 2 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 2 >;

    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 3 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 3 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 3 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 3 >;
    template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3 >;

  } // namespace Fem

} // namespace Dune
