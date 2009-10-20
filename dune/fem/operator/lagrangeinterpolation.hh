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
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of local functions
    typedef typename DiscreteFunctionType::LocalFunctionType
      LocalFunctionType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    //! type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType::LagrangePointSetType
      LagrangePointSetType;
    //! type of vectors in function's domain
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! type of vectors in function's range
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  public:
    /** interpolate an analytical function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  globally) at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class FunctionType >
    static void interpolateFunction ( const FunctionType &function,
                                      DiscreteFunctionType &discreteFunction );

  private:
    template< class FunctionType, bool hasLocalFunction >
    struct CallInterpolateDiscreteFunction;

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
    static void
    interpolateDiscreteFunction ( const FunctionType &function,
                                  DiscreteFunctionType &discreteFunction );
  };



  template< class DiscreteFunctionImp >    
  template< class FunctionType >
  inline void LagrangeInterpolation< DiscreteFunctionImp >
    :: interpolateFunction ( const FunctionType &function,
                             DiscreteFunctionType &discreteFunction )
  {
    const bool hasLocalFunction = Conversion< FunctionType, HasLocalFunction >::exists;
    CallInterpolateDiscreteFunction< FunctionType, hasLocalFunction >::call( function, discreteFunction );
  }



  template< class DiscreteFunctionType >
  template< class FunctionType >
  struct LagrangeInterpolation< DiscreteFunctionType >
    ::CallInterpolateDiscreteFunction< FunctionType, true >
  {
    static void call( const FunctionType &function,
                      DiscreteFunctionType &discreteFunction )
    {
      interpolateDiscreteFunction( function, discreteFunction );
    }
  };

  template< class DiscreteFunctionType >
  template< class FunctionType >
  struct LagrangeInterpolation< DiscreteFunctionType >
    ::CallInterpolateDiscreteFunction< FunctionType, false >
  {
    static void call( const FunctionType &function,
                      DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef DiscreteFunctionAdapter< FunctionType, GridPartType >
        DiscreteFunctionAdapterType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();
      DiscreteFunctionAdapterType dfAdapter( "function", function, dfSpace.gridPart() );
      interpolateDiscreteFunction( dfAdapter, discreteFunction );
    }
  };


  
  template< class DiscreteFunctionImp >
  template< class FunctionType >
  inline void LagrangeInterpolation< DiscreteFunctionImp >
    ::interpolateDiscreteFunction ( const FunctionType &function,
                                    DiscreteFunctionType &discreteFunction )
  {
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    typedef typename FunctionType::LocalFunctionType FunctionLocalFunctionType;

    const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

    IteratorType endit = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != endit; ++it )
    {
      const LagrangePointSetType &lagrangePointSet
        = dfSpace.lagrangePointSet( *it );

      FunctionLocalFunctionType f_local = function.localFunction( *it );
      LocalFunctionType df_local = discreteFunction.localFunction( *it );
      
      // assume point based local dofs 
      const int numDofs = df_local.numDofs();
      for( int i = 0; i < numDofs; ++i )
      {
        RangeType phi;
        f_local.evaluate( lagrangePointSet[ i/dimRange ], phi );
        df_local[ i ] = phi[ i % dimRange ];
      }
    }
  }

}

#endif // #ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
