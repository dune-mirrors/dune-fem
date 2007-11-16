namespace Dune 
{

  template< class DiscreteFunctionImp >    
  template< class FunctionType >
  inline void LagrangeInterpolation< DiscreteFunctionImp >
    :: interpolateFunction ( const FunctionType &function,
                             DiscreteFunctionType &discreteFunction )
  {
    enum { hasLocalFunction = Conversion< FunctionType, HasLocalFunction > :: exists  };
    
    callInterpolateDiscreteFunction< FunctionType, hasLocalFunction >
      :: call( function, discreteFunction );
  }



  template< class DiscreteFunctionType >
  template< class FunctionType >
  struct LagrangeInterpolation< DiscreteFunctionType >
    :: callInterpolateDiscreteFunction< FunctionType, true >
  {
    inline static void call( const FunctionType &function,
                             DiscreteFunctionType &discreteFunction )
    {
      interpolateDiscreteFunction( function, discreteFunction );
    }
  };

  template< class DiscreteFunctionType >
  template< class FunctionType >
  struct LagrangeInterpolation< DiscreteFunctionType >
    :: callInterpolateDiscreteFunction< FunctionType, false >
  {
    inline static void call( const FunctionType &function,
                             DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
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
    :: interpolateDiscreteFunction ( const FunctionType &function,
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



#if DUNE_FEM_COMPATIBILITY
  template< class DiscreteFunctionImp >
  template< class EntityFunctionType >
  inline void LagrangeInterpolation< DiscreteFunctionImp >
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
#endif
  
}
