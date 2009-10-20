#include <config.h>

#include "basefunctions.hh"

namespace Dune
{

  namespace Library
  {

    template< class Field, unsigned int dim, unsigned int pOrder >
    void Instantiate_LagrangeBaseFunctionFactory ()
    {
      LagrangeBaseFunctionFactory
        < FunctionSpace< Field, Field, dim, 1 >, dim, pOrder >
        factory( GeometryType( GeometryType :: simplex, dim ) );
      factory.numBaseFunctions();
    }

    void Instantiate_LagrangeSpace_Templates ()
    {
      Instantiate_LagrangeBaseFunctionFactory< double, 1, 1 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 2, 1 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 3, 1 >();

      Instantiate_LagrangeBaseFunctionFactory< float, 1, 1 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 2, 1 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 3, 1 >();
      
      Instantiate_LagrangeBaseFunctionFactory< double, 1, 2 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 2, 2 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 3, 2 >();
      
      Instantiate_LagrangeBaseFunctionFactory< float, 1, 2 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 2, 2 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 3, 2 >();

      Instantiate_LagrangeBaseFunctionFactory< double, 1, 3 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 2, 3 >();
      Instantiate_LagrangeBaseFunctionFactory< double, 3, 3 >();
      
      Instantiate_LagrangeBaseFunctionFactory< float, 1, 3 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 2, 3 >();
      Instantiate_LagrangeBaseFunctionFactory< float, 3, 3 >();
    }
  }
  
}
