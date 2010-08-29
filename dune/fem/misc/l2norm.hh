#ifndef DUNE_FEM_L2NORM_HH
#define DUNE_FEM_L2NORM_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

namespace Dune
{

  // L2Norm
  // ------

  template< class GridPart >
  class L2Norm
  {
    typedef L2Norm< GridPart > ThisType;

  public:
    typedef GridPart GridPartType;

  protected:
    template< class Function >
    struct FunctionSquare;

    template< class UFunction, class VFunction >
    struct FunctionDistance;

    typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
    typedef typename GridIteratorType::Entity EntityType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef Integrator< QuadratureType > IntegratorType;

  public:
    explicit L2Norm ( const GridPartType &gridPart );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

  protected:
    const GridPartType &gridPart () const { return gridPart_; }

    typename GridPartType::GridType::CollectiveCommunication comm () const
    {
      return gridPart().grid().comm();
    }

  private:
    const GridPartType &gridPart_;
  };



  // WeightedL2Norm
  // --------------
  
  template< class WeightFunction >
  class WeightedL2Norm
  : public L2Norm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType >
  {
    typedef WeightedL2Norm< WeightFunction > ThisType;
    typedef L2Norm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType > BaseType;

  public:
    typedef WeightFunction WeightFunctionType;

    typedef typename WeightFunctionType::DiscreteFunctionSpaceType WeightFunctionSpaceType;
    typedef typename WeightFunctionSpaceType::GridPartType GridPartType;
   
  protected:
    template< class Function >
    struct WeightedFunctionSquare;
    
    typedef typename WeightFunctionType::LocalFunctionType LocalWeightFunctionType;
    typedef typename WeightFunctionType::RangeType WeightType;
    
    typedef typename BaseType::GridIteratorType GridIteratorType;
    typedef typename BaseType::IntegratorType IntegratorType;

    typedef typename GridIteratorType::Entity EntityType;

    using BaseType::gridPart;
    using BaseType::comm;

  public:
    explicit WeightedL2Norm ( const WeightFunctionType &weightFunction );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

  private:
    const WeightFunctionType &weightFunction_;
  };



  // Implementation of L2Norm
  // ------------------------
  
  template< class GridPart >
  inline L2Norm< GridPart >::L2Norm ( const GridPartType &gridPart )
  : gridPart_( gridPart )
  {}


  template< class GridPart >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType::RangeFieldType
  L2Norm< GridPart >::norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

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
  inline typename UDiscreteFunctionType::RangeFieldType
  L2Norm< GridPart >
    ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

    typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = u.space().order();
    const unsigned int vorder = v.space().order();
    const unsigned int order = 2 * std::max( uorder, vorder );
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
  struct L2Norm< GridPart >::FunctionSquare
  {
    typedef Function FunctionType;

    typedef typename FunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

    explicit FunctionSquare ( const FunctionType &function )
    : function_( function )
    {}
    
    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      typename FunctionType::RangeType phi;
      function_.evaluate( x, phi );
      ret = phi * phi;
    }

  private:
    const FunctionType &function_;
  };


  template< class GridPart >
  template< class UFunction, class VFunction >
  struct L2Norm< GridPart >::FunctionDistance
  {
    typedef UFunction UFunctionType;
    typedef VFunction VFunctionType;

    typedef typename UFunctionType::RangeFieldType RangeFieldType;
    typedef typename UFunctionType::RangeType RangeType;
    typedef typename UFunctionType::JacobianRangeType JacobianRangeType;

    FunctionDistance ( const UFunctionType &u, const VFunctionType &v )
    : u_( u ), v_( v )
    {}
    
    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      RangeType phi;
      u_.evaluate( x, ret );
      v_.evaluate( x, phi );
      ret -= phi;
    }

    template< class Point >
    void jacobian ( const Point &x, JacobianRangeType &ret ) const
    {
      JacobianRangeType phi;
      u_.jacobian( x, ret );
      v_.jacobian( x, phi );
      ret -= phi;
    }

  private:
    const UFunctionType &u_;
    const VFunctionType &v_;
  };



  // Implementation of WeightedL2Norm
  // --------------------------------
  
  template< class WeightFunction >
  inline WeightedL2Norm< WeightFunction >
    ::WeightedL2Norm ( const WeightFunctionType &weightFunction )
  : BaseType( weightFunction.space().gridPart ),
    weightFunction_( weightFunction )
  {
    dune_static_assert( (WeightFunctionSpaceType::dimRange == 1),
                        "Wight function must be scalar." );
  }


  template< class WeightFunction >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType::RangeFieldType
  WeightedL2Norm< WeightFunction >
    ::norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

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
  inline typename UDiscreteFunctionType::RangeFieldType
  WeightedL2Norm< WeightFunction >
    ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

    typedef typename BaseType::template FunctionDistance
      < ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = u.space().order();
    const unsigned int vorder = v.space().order();
    const unsigned int order = 2 * std::max( uorder, vorder )
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
  struct WeightedL2Norm< WeightFunction >::WeightedFunctionSquare
  {
    typedef Function FunctionType;

    typedef typename FunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

    WeightedFunctionSquare ( const LocalWeightFunctionType &weightFunction,
                             const FunctionType &function )
    : weightFunction_( weightFunction ),
      function_( function )
    {}
    
    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      WeightType weight;
      weightFunction_.evaluate( x, weight );

      typename FunctionType::RangeType phi;
      function_.evaluate( x, phi );
      ret = weight[ 0 ] * (phi * phi);
    }

  private:
    const LocalWeightFunctionType &weightFunction_;
    const FunctionType &function_;
  };

}

#endif // #ifndef DUNE_FEM_L2NORM_HH
