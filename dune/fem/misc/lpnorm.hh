#ifndef DUNE_FEM_LPNORM_HH
#define DUNE_FEM_LPNORM_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

namespace Dune
{

  // TODO weighte LP norm might be adapted later
  // LPNorm
  //
  // !!!!! It is not cleared which quadratur order have to be applied for p > 2!!!!
  // !!!!! For p = 1 this norm does not work !!!!
  // ------


  //! Quadrature Order Interface
  struct OrderCalculatorInterface
  {
    virtual int operator() (const double p)=0;
  };

  //! default Quadrature Order Calculator
  //  can be re-implemented in order to use 
  //  a different type of calculation
  //  which can be sepcified in the second template argument of LPNorm
  struct DefaultOrderCalculator : public OrderCalculatorInterface
  {
    int operator() (const double p)
    {
      int ret=0;
      const double q = p / (p-1);
      double max = std::max(p,q);
      assert(max < std::numeric_limits<int>::max()/2. );
      ret = max +1;
      return ret;
    }
  };

  template< class GridPart, class OrderCalculator = DefaultOrderCalculator >
  class LPNorm
  {
    typedef LPNorm< GridPart > ThisType;

  public:
    typedef GridPart GridPartType;

  protected:
    template< class Function >
    struct FunctionMultiplicator;

    template< class UFunction, class VFunction >
    struct FunctionDistance;

    typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
    typedef typename GridIteratorType::Entity EntityType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef Integrator< QuadratureType > IntegratorType;

  public:
    explicit LPNorm ( const GridPartType &gridPart, const double p );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

  protected:
    const GridPartType &gridPart () const { return gridPart_; }

    typename GridPartType::GridType::Traits::CollectiveCommunication comm () const
    {
      return gridPart().grid().comm();
    }

    int order ( const int spaceOrder ) const ;

  private:
    const GridPartType &gridPart_;
    double p_ ;
    OrderCalculator *orderCalc;
  };


#if 0
  // WeightedLPNorm
  // --------------
  
  template< class WeightFunction >
  class WeightedLPNorm
  : public LPNorm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType >
  {
    typedef WeightedLPNorm< WeightFunction > ThisType;
    typedef LPNorm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType > BaseType;

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
    explicit WeightedLPNorm ( const WeightFunctionType &weightFunction, const double p_ );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

  private:
    const WeightFunctionType &weightFunction_;
    const double p_;
  };

#endif

  // Implementation of LPNorm
  // ------------------------
 
  template< class GridPart, class OrderCalculator >
  inline LPNorm< GridPart, OrderCalculator >::LPNorm ( const GridPartType &gridPart, const double p )
  : gridPart_( gridPart ),
    p_(p)
  {
    orderCalc = new OrderCalculator();
  }

  template< class GridPart, class OrderCalculator>
  inline int LPNorm< GridPart, OrderCalculator>::order(const int spaceOrder) const
  {
    return spaceOrder * orderCalc->operator() (p_);
  }


  template< class GridPart, class OrderCalculator>
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType::RangeFieldType
  LPNorm< GridPart, OrderCalculator >::norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    unsigned int order = order( u.space() );
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionMultiplicator< LocalFunctionType > ulocalp( ulocal, p_ );

      integrator.integrateAdd( entity, ulocalp, sum );
    }

    return std::pow (  comm().sum( sum[ 0 ] ), 1./p_ );
  }

  template< class GridPart, class OrderCalculator >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType::RangeFieldType
  LPNorm< GridPart, OrderCalculator >
    ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

    typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = order( u.space().order() );
    const unsigned int vorder = order( v.space().order() );
    const unsigned int order = std::max( uorder, vorder );
    IntegratorType integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionMultiplicator< LocalDistanceType > distp( dist, p_ );
      
      integrator.integrateAdd( entity, distp, sum );
    }

    return std::pow( comm().sum( sum[ 0 ] ), 1./p_ );
  }

 
  template< class GridPart, class OrderCalculator >
  template< class Function >
  struct LPNorm< GridPart, OrderCalculator >::FunctionMultiplicator
  {
    typedef Function FunctionType;

    typedef typename FunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

    explicit FunctionMultiplicator ( const FunctionType &function, double p )
    : function_( function ),
      p_(p)
    {}
    
    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      typename FunctionType::RangeType phi;
      function_.evaluate( x, phi );
      ret = std :: pow ( phi.two_norm(), p_);
    }

  private:
    const FunctionType &function_;
    double p_;
  };


  template< class GridPart, class OrderCalculator >    
  template< class UFunction, class VFunction >
  struct LPNorm< GridPart, OrderCalculator >::FunctionDistance
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


#if 0
  // Implementation of WeightedL2Norm
  // --------------------------------
  
  template< class WeightFunction >
  inline WeightedLPNorm< WeightFunction >
    ::WeightedLPNorm ( const WeightFunctionType &weightFunction )
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
#endif
}

#endif // #ifndef DUNE_FEM_L2NORM_HH
