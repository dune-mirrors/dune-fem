#ifndef DUNE_LAGRANGEDOFHANDLER_HH
#define DUNE_LAGRANGEDOFHANDLER_HH

#warning "dune/fem/operator/lagrangedofhandler.hh is deprecated."

#include <dune/fem/space/lagrangespace/dofhandler.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

namespace Dune 
{

  /** \class LagrangeInterpolator
   *  \brief Generates the Lagrange Interpolation of an analytical function
   */
  class LagrangeInterpolator
  {
  public:
    /** interpolate an analytical function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  globally) at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \note To interpolate discrete functions, use interpolateEntityFunction.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class FunctionType, class DiscreteFunctionType >
    void interpolFunction ( FunctionType &function,
                            DiscreteFunctionType &discreteFunction )
    {
      LagrangeInterpolation< DiscreteFunctionType >
      :: interpolateFunction( function, discreteFunction );
    }

    /** interpolate an entity-function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  locally at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class EntityFunctionType, class DiscreteFunctionType >
    static void interpolEntityFunction ( EntityFunctionType &function,
                                         DiscreteFunctionType &discreteFunction )
    {
      LagrangeInterpolation< DiscreteFunctionType >
      :: interpolateEntityFunction( function, discreteFunction );
    }
  };

}

#endif
