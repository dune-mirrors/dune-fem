namespace Dune
{

  // L2 Norm
  // -------
  
  template< class GridPart >
  inline L2Norm< GridPart > :: L2Norm ( const GridPartType &gridPart )
  : gridPart_( gridPart )
  {
  }



  template< class GridPart >
  inline L2Norm< GridPart > :: L2Norm ( const ThisType &other )
  : gridPart_( other.gridPart_ )
  {
  }


  
  template< class GridPart >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType :: RangeFieldType
  L2Norm< GridPart > :: norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionType :: RangeType RangeType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    RangeFieldType sum( 0 );
    unsigned int order = 2 * u.space().order();

    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;
      //const GeometryType &geometry = entity.geometry();

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionSquare< LocalFunctionType > ulocal2( ulocal );

      FieldVector< RangeFieldType, 1 > localSum;
      QuadratureType quadrature( entity, order );
      quadrature.integrate( ulocal2, localSum );

      sum += localSum[ 0 ];
    }

    return sqrt( comm().sum( sum ) );
  }

 
  
  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType :: RangeFieldType
  L2Norm< GridPart > :: distance ( const UDiscreteFunctionType &u,
                                   const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType :: RangeFieldType RangeFieldType;
    typedef typename UDiscreteFunctionType :: RangeType RangeType;

    typedef typename UDiscreteFunctionType :: LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType :: LocalFunctionType VLocalFunctionType;

    typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    RangeFieldType sum( 0 );
    const unsigned int uorder = 2 * u.space().order();
    const unsigned int vorder = 2 * v.space().order();
    const unsigned int order = std :: max( uorder, vorder );

    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionSquare< LocalDistanceType > dist2( dist );
      
      FieldVector< RangeFieldType, 1 > localSum;
      QuadratureType quadrature( entity, order );
      quadrature.integrate( dist2, localSum );

      sum += localSum[ 0 ];
    }

    return sqrt( comm().sum( sum ) );
  }


  
  template< class GridPart >
  template< class Function >
  class L2Norm< GridPart > :: FunctionSquare
  {
  public:
    typedef Function FunctionType;

    typedef typename FunctionType :: RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

  protected:
    const FunctionType &function_;
  
  public:
    inline explicit FunctionSquare ( const FunctionType &function )
    : function_( function )
    {}
    
    template< class Point >
    inline void evaluate ( const Point &x,
                           RangeType &ret ) const
    {
      typename FunctionType :: RangeType phi;
      function_.evaluate( x, phi );
      ret = phi * phi;
    }
  };



  template< class GridPart >
  template< class UFunction, class VFunction >
  class L2Norm< GridPart > :: FunctionDistance
  {
  public:
    typedef UFunction UFunctionType;
    typedef VFunction VFunctionType;

    typedef typename UFunctionType :: RangeFieldType RangeFieldType;
    typedef typename UFunctionType :: RangeType RangeType;

  protected:
    const UFunctionType &u_;
    const VFunctionType &v_;
  
  public:
    inline explicit FunctionDistance ( const UFunctionType &u,
                                       const VFunctionType &v )
    : u_( u ),
      v_( v )
    {}
    
    template< class Point >
    inline void evaluate ( const Point &x,
                           RangeType &ret ) const
    {
      RangeType phi;
      u_.evaluate( x, ret );
      v_.evaluate( x, phi );
      ret -= phi;
    }
  };



  // Weighted L2 Norm
  // ----------------
  
  template< class WeightFunction >
  inline WeightedL2Norm< WeightFunction >
    :: WeightedL2Norm ( const WeightFunctionType &weightFunction )
  : BaseType( weightFunction.space().gridPart ),
    weightFunction_( weightFunction )
  {
    CompileTimeChecker< (WeightFunctionSpaceType :: dimRange == 1) >
      __WEIGHTFUNCTION_MUST_BE_SCALAR__;
  }



  template< class WeightFunction >
  inline WeightedL2Norm< WeightFunction >
    :: WeightedL2Norm ( const ThisType &other )
  : BaseType( other ),
    weightFunction_( other.weightFunction_ )
  {
  }


  
  template< class WeightFunction >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType :: RangeFieldType
  WeightedL2Norm< WeightFunction >
    :: norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionType :: RangeType RangeType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    RangeFieldType sum( 0 );
    unsigned int order = 2 * u.space().order() + weightFunction_.space().order();

    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      LocalFunctionType ulocal = u.localFunction( entity );
     
      WeightedFunctionSquare< LocalFunctionType > ulocal2( wflocal, ulocal );

      FieldVector< RangeFieldType, 1 > localSum;
      QuadratureType quadrature( entity, order );
      quadrature.integrate( ulocal2, localSum );

      sum += localSum[ 0 ];
    }

    return sqrt( comm().sum( sum ) );
  }

 
  
  template< class WeightFunction >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType :: RangeFieldType
  WeightedL2Norm< WeightFunction >
    :: distance ( const UDiscreteFunctionType &u,
                  const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType :: RangeFieldType RangeFieldType;
    typedef typename UDiscreteFunctionType :: RangeType RangeType;

    typedef typename UDiscreteFunctionType :: LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType :: LocalFunctionType VLocalFunctionType;

    typedef typename BaseType :: template FunctionDistance
      < ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    RangeFieldType sum( 0 );
    const unsigned int uorder = 2 * u.space().order();
    const unsigned int vorder = 2 * v.space().order();
    const unsigned int order = std :: max( uorder, vorder )
                               + weightFunction_.space().order();

    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );
     
      LocalDistanceType dist( ulocal, vlocal );
      WeightedFunctionSquare< LocalDistanceType > dist2( wflocal, dist );
      
      FieldVector< RangeFieldType, 1 > localSum;
      QuadratureType quadrature( entity, order );
      quadrature.integrate( dist2, localSum );

      sum += localSum[ 0 ];
    }

    return sqrt( comm().sum( sum ) );
  }


  
  template< class WeightFunction >
  template< class Function >
  class WeightedL2Norm< WeightFunction > :: WeightedFunctionSquare
  {
  public:
    typedef Function FunctionType;

    typedef typename FunctionType :: RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

  protected:
    const LocalWeightFunctionType &weightFunction_;
    const FunctionType &function_;
  
  public:
    inline WeightedFunctionSquare ( const LocalWeightFunctionType &weightFunction,
                                    const FunctionType &function )
    : weightFunction_( weightFunction ),
      function_( function )
    {}
    
    template< class Point >
    inline void evaluate ( const Point &x,
                           RangeType &ret ) const
    {
      WeightType weight;
      weightFunction_.evaluate( x, weight );

      typename FunctionType :: RangeType phi;
      function_.evaluate( x, phi );
      ret = weight[ 0 ] * (phi * phi);
    }
  };


}
