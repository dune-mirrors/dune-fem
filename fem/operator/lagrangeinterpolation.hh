#ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
#define DUNE_FEM_LAGRANGEINTERPOLATION_HH

#include <dune/common/typetraits.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>

namespace Dune 
{

  /** \class LagrangeInterpolation
   *  \brief Generates the Lagrange Interpolation of an analytic function
   */
  template< class DiscreteFunctionImp >
  class LagrangeInterpolation
  {
  public:
    //! type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of local functions
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
      LagrangePointSetType;
    //! type of vectors in function's domain
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of vectors in function's range
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
    static void interpolateFunction ( const FunctionType &function,
                                      DiscreteFunctionType &discreteFunction )
    {
      enum { hasLocalFunction = Conversion< FunctionType, HasLocalFunction > :: exists > };
      
      callInterpolateDiscreteFunction< FunctionType, hasLocalFunction >
        :: call( function, discreteFunction );
    }

  private:
    template< class FunctionType, bool hasLocalFunction >
    struct callInterpolateDiscreteFunction;

    template< class FunctionType >
    struct callInterpolateDiscreteFunction< FunctionType, true >
    {
      static void call( const FunctionType &function, 
                        DiscreteFunctionType &discreteFunction )
      {
        interpolateDiscreteFunction( function, discreteFunction );
      }
    };

    template< class FunctionType >
    struct callInterpolateDiscreteFunction< FunctionType, false >
    {
      static void call( const FunctionType &function,
                        DiscreteFunctionType &discreteFunction )
      {
        typedef DiscreteFunctionAdapter< FunctionType, GridPartType >
          DiscreteFunctionAdapterType;

        const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

        DiscreteFunctionAdapterType dfAdapter( "function", function, dfSpace.gridPart() );

        interpolateDiscreteFunction( dfAdapter, discreteFunction );
      }
    };

  protected:
    /** interpolate a discrete function into a Lagrange discrete function
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
    template< class FunctionType >
    static void interpolateDiscreteFunction
                ( const FunctionType &function,
                  DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      typedef typename FunctionType :: LocalFunctionType FunctionLocalFunctionType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

      //discreteFunction.clear();

      IteratorType endit = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != endit; ++it )
      {
        const LagrangePointSetType &lagrangePointSet
          = dfSpace.lagrangePointSet( *it );

        FunctionLocalFunctionType f_local = function.localFunction( *it );
        LocalFunctionType df_local = discreteFunction.localFunction( *it );
        
        const int numDofs = df_local.numDofs();
        for( int i = 0; i < numDofs; ++i )
        {
          RangeType phi;
          f_local.evaluate( lagrangePointSet, i, phi );
          df_local[ i ] = phi[ 0 ];
        }
      }
    }


  public:
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
    static void interpolateEntityFunction ( EntityFunctionType &function,
                                            DiscreteFunctionType &discreteFunction ) DUNE_DEPRECATED;
  };


 
  template< class DiscreteFunctionImp >
  template< class EntityFunctionType >
  void LagrangeInterpolation< DiscreteFunctionImp >
    :: interpolateEntityFunction ( EntityFunctionType &function,
                                   DiscreteFunctionType &discreteFunction )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    const DiscreteFunctionSpaceType &functionSpace
      = discreteFunction.space();

    //discreteFunction.clear();

    IteratorType endit = functionSpace.end();
    for( IteratorType it = functionSpace.begin(); it != endit; ++it )
    {
      const LagrangePointSetType &lagrangePointSet
        = functionSpace.lagrangePointSet( *it );

      LocalFunctionType localFunction = discreteFunction.localFunction( *it );
      function.init( *it );
      
      const int numDofs = localFunction.numDofs();
      for( int i = 0; i < numDofs; ++i )
      {
        RangeType phi;
        function.evaluate( lagrangePointSet.point( i ), phi );
        localFunction[ i ] = phi[ 0 ];
      }
    }
  }
  
}

#endif
