#ifndef DUNE_FEM_L2NORM_HH
#define DUNE_FEM_L2NORM_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

namespace Dune
{

  template< class GridPart >
  class L2Norm
  {
  public:
    typedef GridPart GridPartType;

  private:
    typedef L2Norm< GridPartType > ThisType;

  protected:
    template< class Function >
    class FunctionSquare;

    template< class UFunction, class VFunction >
    class FunctionDistance;

  protected:
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      GridIteratorType;
    typedef typename GridIteratorType :: Entity EntityType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef Integrator< QuadratureType > IntegratorType;

  protected:
    const GridPartType &gridPart_;

  public:
    inline explicit L2Norm ( const GridPartType &gridPart );
    inline L2Norm ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );

  public:
    template< class DiscreteFunctionType >
    inline typename DiscreteFunctionType :: RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename UDiscreteFunctionType :: RangeFieldType
    distance ( const UDiscreteFunctionType &u,
               const VDiscreteFunctionType &v ) const;

  protected:
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

    typename GridPartType::GridType ::CollectiveCommunication comm () const
    {
      return gridPart().grid().comm();
    }
  };


  
  template< class WeightFunction >
  class WeightedL2Norm
  : public L2Norm
    < typename WeightFunction :: DiscreteFunctionSpaceType :: GridPartType >
  {
  public:
    typedef WeightFunction WeightFunctionType;

   public:
    typedef typename WeightFunctionType :: DiscreteFunctionSpaceType
      WeightFunctionSpaceType;
    typedef typename WeightFunctionSpaceType :: GridPartType GridPartType;
   
  private:
    typedef WeightedL2Norm< WeightFunctionType > ThisType;
    typedef L2Norm< GridPartType > BaseType;

  protected:
    template< class Function >
    class WeightedFunctionSquare;
    
  protected:
    typedef typename WeightFunctionType :: LocalFunctionType
      LocalWeightFunctionType;
    typedef typename WeightFunctionType :: RangeType WeightType;
    
    typedef typename BaseType :: GridIteratorType GridIteratorType;
    typedef typename BaseType :: IntegratorType IntegratorType;

    typedef typename GridIteratorType :: Entity EntityType;

  protected:
    const WeightFunctionType &weightFunction_;

  protected:
    using BaseType :: gridPart;
    using BaseType :: comm;

  public:
    inline explicit WeightedL2Norm ( const WeightFunctionType &weightFunction );
    inline WeightedL2Norm ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );

  public:
    template< class DiscreteFunctionType >
    inline typename DiscreteFunctionType :: RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename UDiscreteFunctionType :: RangeFieldType
    distance ( const UDiscreteFunctionType &u,
               const VDiscreteFunctionType &v ) const;
  };



  // Implementation of L2Norm
  // ------------------------
  
  template< class GridPart >
  inline L2Norm< GridPart > :: L2Norm ( const GridPartType &gridPart )
  : gridPart_( gridPart )
  {}


  template< class GridPart >
  inline L2Norm< GridPart > :: L2Norm ( const ThisType &other )
  : gridPart_( other.gridPart_ )
  {}

  
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



  // Implementation of WeightedL2Norm
  // --------------------------------
  
  template< class WeightFunction >
  inline WeightedL2Norm< WeightFunction >
    :: WeightedL2Norm ( const WeightFunctionType &weightFunction )
  : BaseType( weightFunction.space().gridPart ),
    weightFunction_( weightFunction )
  {
    dune_static_assert( (WeightFunctionSpaceType::dimRange == 1),
                        "Wight function must be scalar." );
  }


  template< class WeightFunction >
  inline WeightedL2Norm< WeightFunction >
    :: WeightedL2Norm ( const ThisType &other )
  : BaseType( other ),
    weightFunction_( other.weightFunction_ )
  {}

  
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

#endif // #ifndef DUNE_FEM_L2NORM_HH
