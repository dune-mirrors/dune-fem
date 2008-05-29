#ifndef DUNE_FEM_SOURCEPROJECTION_HH
#define DUNE_FEM_SOURCEPROJECTION_HH

#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/operator/integrationoperator.hh>

namespace Dune
{

  template< class DomainFunctionImp, class RangeFunctionImp >
  class DefaultEllipticSourceProjectionTraits
  {
  public:
    typedef DomainFunctionImp DomainFunctionType;
    typedef RangeFunctionImp RangeFunctionType;

  private:
    typedef DefaultEllipticSourceProjectionTraits< DomainFunctionImp, RangeFunctionImp >
      ThisType;

  public:
    typedef typename RangeFunctionImp :: FunctionSpaceType DiscreteFunctionSpaceType;
   
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
  };



  template< class TraitsImp >
  class LocalEllipticSourceProjection
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef LocalEllipticSourceProjection< TraitsType > ThisType;

  public:
    typedef typename TraitsType :: DomainFunctionType DomainFunctionType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename DomainFunctionType :: RangeFieldType DomainFieldType;
    typedef typename RangeFunctionType :: RangeFieldType RangeFieldType;

    typedef typename DomainFunctionType :: FunctionSpaceType DomainFunctionSpaceType;
    typedef typename RangeFunctionType :: FunctionSpaceType RangeFunctionSpaceType;

    typedef typename TraitsType :: QuadratureType QuadratureType;
    typedef typename QuadratureType :: CoordinateType QuadraturePointType;

    typedef typename TraitsType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

  public:
    inline LocalEllipticSourceProjection ()
    {
    }

    inline LocalEllipticSourceProjection ( const ThisType &other )
    {
    }
    
    template< class EntityType, class T >
    inline void operator() ( const EntityType &entity,
                             const DomainFunctionType &u,
                             LocalFunction< T > &w ) const
    {
      typedef LocalFunction< T > RangeLocalFunctionType;

      // geometry type for the entity
      typedef typename EntityType :: Geometry GeometryType;

      // type of local function for domain function
      typedef typename DomainFunctionType :: LocalFunctionType
        DomainLocalFunctionType;
      typedef typename DomainFunctionType :: RangeType DomainFunctionRangeType;

      // type of base function set
      typedef typename RangeLocalFunctionType :: BaseFunctionSetType
        BaseFunctionSetType;
      typedef typename RangeLocalFunctionType :: RangeType RangeFunctionRangeType;

      // obtain geometry from the entity
      const GeometryType &geometry = entity.geometry();
      
      // obtain a local function for the argument function
      DomainLocalFunctionType u_local = u.localFunction( entity );

      // obtain range function's base function set
      const BaseFunctionSetType &baseFunctionSet = w.baseFunctionSet();
      const unsigned int numDofs = w.numDofs();

      // clear the destination function
      for( unsigned int i = 0; i < numDofs; ++i )
        w[ i ] = 0;

      // loop adding up the destination function
      QuadratureType quadrature( entity, 2 * polynomialOrder + 2 );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // get quadrature point
        const QuadraturePointType &point = quadrature.point( pt );

        // weight of this point in the integral
        const RangeFieldType weight
          = geometry.integrationElement( point ) * quadrature.weight( pt );

        // evaluate argument function
        DomainFunctionRangeType phi;
        u_local.evaluate( quadrature[ pt ], phi );
        
        // multiply with all base functions
        for( unsigned int i = 0; i < numDofs; ++i )
        {
          RangeFunctionRangeType psi;
          baseFunctionSet.evaluate( i, quadrature[ pt ], psi );

          // update destination function
          w[ i ] += weight * (phi * psi);
        }
      }
    }
  };



  template< class TraitsImp >
  class EllipticSourceProjection
  : public IntegrationOperator
    < DefaultIntegrationOperatorTraits
      < LocalEllipticSourceProjection< TraitsImp >,
        typename TraitsImp :: RangeFunctionType
      >,
      false
    >
  {
  public:
    typedef TraitsImp TraitsType;

    typedef LocalEllipticSourceProjection< TraitsType > LocalOperatorType;

    typedef typename TraitsType :: DomainFunctionType DomainFunctionType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

  private:
    typedef DefaultIntegrationOperatorTraits< LocalOperatorType, RangeFunctionType >
      IntegrationOperatorTraitsType;

    typedef EllipticSourceProjection< TraitsType > ThisType;
    typedef IntegrationOperator< IntegrationOperatorTraitsType, false > BaseType;

  protected:
    LocalOperatorType localOperator_;

  public:
    inline EllipticSourceProjection ()
    : BaseType( localOperator_ ),
      localOperator_()
    {
    }
  };

}

#endif
