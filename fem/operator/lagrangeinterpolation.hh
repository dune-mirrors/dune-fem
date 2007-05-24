#ifndef DUNE_LAGRANGEINTERPOLATION_HH
#define DUNE_LAGRANGEINTERPOLATION_HH

#include <dune/fem/misc/entityfunction.hh>

namespace Dune 
{

  /** \class LagrangeInterpolation
   *  \brief Generates the Lagrange Interpolation of an analytic function
   */
  template< class DiscreteFunctionImp >
  class LagrangeInterpolation
  {
  public:
    //! Type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! Type of discrete function space
    typedef typename DiscreteFunctionType :: FunctionSpaceType
      DiscreteFunctionSpaceType;
    //! Type of local functions
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;
    //! Type of grid
    typedef typename DiscreteFunctionType :: GridType GridType;

    //! Type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
      LagrangePointSetType;
    //! Type of vectors in function's domain
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! Type of vectors in function's range
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

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
    template< class FunctionType >
    static void interpolateFunction ( FunctionType &function,
                                      DiscreteFunctionType &discreteFunction )
    {
      typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;

      EntityFunctionAdapter< Entity0Type, FunctionType >
        entityFunction( function );

      interpolateEntityFunction( entityFunction, discreteFunction );
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
    template< class EntityFunctionType >
    static void interpolateEntityFunction
                ( EntityFunctionType &function,
                  DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      const DiscreteFunctionSpaceType &functionSpace
        = discreteFunction.space();

      discreteFunction.clear();

      IteratorType endit = functionSpace.end();
      for( IteratorType it = functionSpace.begin(); it != endit; ++it ) {
        const LagrangePointSetType &lagrangePointSet
          = functionSpace.lagrangePointSet( *it );

        LocalFunctionType localFunction = discreteFunction.localFunction( *it );
        function.init( *it );
        
        const int numDofs = localFunction.numDofs();
        for( int i = 0; i < numDofs; ++i ) {
          RangeType phi;
          function.evaluate( lagrangePointSet[ i ], phi );
          localFunction[ i ] = phi[ 0 ];
        }
      }
    }
  };

}

#endif
