#ifndef DUNE_FEM_L2NORM_INLINE_HH
#define DUNE_FEM_L2NORM_INLINE_HH

#include "l2norm.hh"

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

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    unsigned int order = 2 * u.space().order();
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionSquare< LocalFunctionType > ulocal2( ulocal );

      integrator.integrateAdd( entity, ulocal2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
  }

 
  
  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType :: RangeFieldType
  L2Norm< GridPart > :: distance ( const UDiscreteFunctionType &u,
                                   const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType :: RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType :: LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType :: LocalFunctionType VLocalFunctionType;

    typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = u.space().order();
    const unsigned int vorder = v.space().order();
    const unsigned int order = 2 * std :: max( uorder, vorder );
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionSquare< LocalDistanceType > dist2( dist );
      
      integrator.integrateAdd( entity, dist2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
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
    typedef typename UFunctionType :: JacobianRangeType JacobianRangeType;

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

    template< class Point >
    inline void jacobian ( const Point &x,
                           JacobianRangeType &ret ) const
    {
      JacobianRangeType phi;
      u_.jacobian( x, ret );
      v_.jacobian( x, phi );
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

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

    unsigned int order = 2 * u.space().order() + weightFunction_.space().order();
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      LocalFunctionType ulocal = u.localFunction( entity );
     
      WeightedFunctionSquare< LocalFunctionType > ulocal2( wflocal, ulocal );

      integrator.integrateAdd( entity, ulocal2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
  }

 
  
  template< class WeightFunction >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType :: RangeFieldType
  WeightedL2Norm< WeightFunction >
    :: distance ( const UDiscreteFunctionType &u,
                  const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType :: RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType :: LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType :: LocalFunctionType VLocalFunctionType;

    typedef typename BaseType :: template FunctionDistance
      < ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = u.space().order();
    const unsigned int vorder = v.space().order();
    const unsigned int order = 2 * std :: max( uorder, vorder )
                               + weightFunction_.space().order();
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );
     
      LocalDistanceType dist( ulocal, vlocal );
      WeightedFunctionSquare< LocalDistanceType > dist2( wflocal, dist );
      
      integrator.integrateAdd( entity, dist2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
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

#endif
