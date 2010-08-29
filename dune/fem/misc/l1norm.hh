#ifndef DUNE_FEM_L1NORM_HH
#define DUNE_FEM_L1NORM_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

namespace Dune
{

  // L1Norm
  // ------

  template< class GridPart >
  class L1Norm
  {
    typedef L1Norm< GridPart > ThisType;

  public:
    typedef GridPart GridPartType;

  protected:
    template< class Function >
    struct FunctionAbs;

    template< class UFunction, class VFunction >
    struct FunctionDistance;

    typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
    typedef typename GridIteratorType::Entity EntityType;
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

  public:
    explicit L1Norm ( const GridPartType &gridPart );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType norm ( const DiscreteFunctionType &u ) const;
    
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



  // Implementation of L1Norm
  // ------------------------
  
  template< class GridPart >
  inline L1Norm< GridPart >::L1Norm ( const GridPartType &gridPart )
  : gridPart_( gridPart )
  {}


  template< class GridPart >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType::RangeFieldType
  L1Norm< GridPart >::norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    unsigned int order = u.space().order();
    Integrator< QuadratureType > integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionAbs< LocalFunctionType > ulocalAbs( ulocal );

      integrator.integrateAdd( entity, ulocalAbs, sum );
    }

    return comm().sum( sum[ 0 ] );
  }

  
  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType::RangeFieldType
  L1Norm< GridPart >
    ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;

    typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

    typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    const unsigned int uorder = u.space().order();
    const unsigned int vorder = v.space().order();
    const unsigned int order = std::max( uorder, vorder );
    Integrator< QuadratureType > integrator( order );

    FieldVector< RangeFieldType, 1 > sum( 0 );
    const GridIteratorType end = gridPart().template end< 0 >();
    for( GridIteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
    {
      const EntityType &entity = *it;

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionAbs< LocalDistanceType > distAbs( dist );
      
      integrator.integrateAdd( entity, distAbs, sum );
    }

    return comm().sum( sum[ 0 ] );
  }

  
  template< class GridPart >
  template< class Function >
  struct L1Norm< GridPart >::FunctionAbs
  {
    typedef Function FunctionType;

    typedef typename FunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

    explicit FunctionAbs ( const FunctionType &function )
    : function_( function )
    {}
    
    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      typename FunctionType::RangeType phi;
      function_.evaluate( x, phi );
      ret = phi.one_norm();
    }

  private:
    const FunctionType &function_;
  };


  template< class GridPart >
  template< class UFunction, class VFunction >
  struct L1Norm< GridPart >::FunctionDistance
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

}

#endif // #ifndef DUNE_FEM_L1NORM_HH
