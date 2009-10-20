#ifndef DUNE_FEM_H1NORM_INLINE_HH
#define DUNE_FEM_H1NORM_INLINE_HH

#include "h1norm.hh"

namespace Dune
{

  // H1 Norm
  // -------
  
  template< class GridPart >
  inline H1Norm< GridPart > :: H1Norm ( const GridPartType &gridPart )
  : BaseType( gridPart )
  {
  }



  template< class GridPart >
  inline H1Norm< GridPart > :: H1Norm ( const ThisType &other )
  : BaseType( other )
  {
  }


  
  template< class GridPart >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType :: RangeFieldType
  H1Norm< GridPart > :: norm ( const DiscreteFunctionType &u ) const
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
      FunctionJacobianSquare< LocalFunctionType > ulocal2( ulocal );

      integrator.integrateAdd( entity, ulocal2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
  }

 
  
  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType :: RangeFieldType
  H1Norm< GridPart > :: distance ( const UDiscreteFunctionType &u,
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
      FunctionJacobianSquare< LocalDistanceType > dist2( dist );

      integrator.integrateAdd( entity, dist2, sum );
    }

    return sqrt( comm().sum( sum[ 0 ] ) );
  }



  template< class GridPart >
  template< class Function >
  class H1Norm< GridPart > :: FunctionJacobianSquare
  {
  public:
    typedef Function FunctionType;

    typedef typename FunctionType :: RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

  protected:
    const FunctionType &function_;
  
  public:
    inline explicit FunctionJacobianSquare ( const FunctionType &function )
    : function_( function )
    {}
    
    template< class Point >
    inline void evaluate ( const Point &x,
                           RangeType &ret ) const
    {
      enum { dimRange = FunctionType :: RangeType :: size };
      
      typename FunctionType :: RangeType phi;
      function_.evaluate( x, phi );
      ret[ 0 ] = phi * phi;

      typename FunctionType :: JacobianRangeType grad;
      function_.jacobian( x, grad );
      for( unsigned int i = 0; i < dimRange; ++i )
        ret[ 0 ] += (grad[ i ] * grad[ i ]);
    }
  };

}

#endif
