#ifndef DUNE_FEM_DIFFUSIONOPERATOR_HH
#define DUNE_FEM_DIFFUSIONOPERATOR_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

//#include <dune/fem/operator/matrix/localmatrix.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/identityoperator.hh>
#include <dune/fem/operator/integrationoperator.hh>
#include <dune/fem/operator/sourceprojection.hh>

#include <dune/fem/operator/model/diffusionmodel.hh>

namespace Dune
{

  template< class SourceFunction, unsigned int polOrder >
  class DefaultDiffusionOperatorTraits
  {
  public:
    typedef SourceFunction SourceFunctionType;

    enum { polynomialOrder = polOrder };

  private:
    typedef DefaultDiffusionOperatorTraits< SourceFunctionType, polynomialOrder >
      ThisType;
    
  public:
    typedef typename SourceFunctionType :: FunctionSpaceType FunctionSpaceType;

    typedef typename SourceFunctionType :: GridPartType GridPartType;

    typedef LagrangeDiscreteFunctionSpace
      < FunctionSpaceType, GridPartType, polynomialOrder, CachingStorage >
      DiscreteFunctionSpaceType;

    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  };



  template< class T, class Model >
  class LocalDiffusionOperator
  {
  public:
    typedef T Traits;

    typedef Model ModelType;

    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

    typedef DiffusionModelInterface< FunctionSpaceType, ModelType >
      ModelInterfaceType;
   
  private:
    typedef LocalDiffusionOperator< Traits, ModelType > ThisType;
 
  public:
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType DomainFunctionSpaceType;
    typedef DiscreteFunctionSpaceType RangeFunctionSpaceType;

    typedef typename Traits :: DiscreteFunctionType DiscreteFunctionType;
    typedef DiscreteFunctionType DomainFunctionType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef typename QuadratureType :: CoordinateType QuadraturePointType;
 
    typedef typename DomainFunctionType :: RangeFieldType DomainFieldType;
    typedef DomainFieldType RangeFieldType;

    enum { polynomialOrder = Traits :: polynomialOrder };

  protected:
    typedef typename EntityType :: Geometry GeometryType;

    typedef FieldMatrix< typename GeometryType :: ctype,
                         GeometryType :: mydimension,
                         GeometryType :: mydimension >
      GeometryJacobianInverseType;

  protected:
    const ModelType &model_;

  public:
    inline explicit LocalDiffusionOperator ( const ModelType &model )
    : model_( model )
    {
    }

    inline LocalDiffusionOperator ( const ThisType &other )
    : model_( other.model_ )
    {
    }
 
    template< class RangeLocalFunction >
    inline void operator() ( const EntityType &entity,
                             const DomainFunctionType &u,
                             RangeLocalFunction &w ) const
    {
      typedef RangeLocalFunction RangeLocalFunctionType;
       // local function type for domain function
      typedef typename DomainFunctionType :: LocalFunctionType
        DomainLocalFunctionType;
    
      // type of base function set
      typedef typename RangeLocalFunctionType :: BaseFunctionSetType
        BaseFunctionSetType;

      // type of jacobians
      typedef typename RangeLocalFunctionType :: RangeFieldType RangeFieldType;
      typedef typename RangeLocalFunctionType :: JacobianRangeType RangeJacobianType;
      typedef typename DomainFunctionSpaceType :: JacobianRangeType DomainJacobianType;

      // obtain geometry from the entity
      const GeometryType &geometry = entity.geometry();

      // obtain local function for argument function
      DomainLocalFunctionType u_local = u.localFunction( entity );

      // obtain base function set
      const BaseFunctionSetType &baseFunctionSet = w.baseFunctionSet();
      const unsigned int numDofs = w.numDofs();

      // clear the destination function
      for( unsigned int i = 0; i < numDofs; ++i )
        w[ i ] = 0;

      // loop adding up the destination function
      QuadratureType quadrature( entity, 2 * polynomialOrder );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // get quadrature point
        const QuadraturePointType &point = quadrature.point( pt );

        // get jacobian inverse of reference mapping
        const GeometryJacobianInverseType &inv
          = geometry.jacobianInverseTransposed( point );

        // weight of this point in the integral
        const RangeFieldType weight
          = geometry.integrationElement( point ) * quadrature.weight( pt );
     
        // evaluate gradient of argument function
        DomainJacobianType u_x;
        u_local.jacobian( quadrature[ pt ], u_x );

        // calculate diffusive flux
        RangeJacobianType diffusiveFlux;
        model_.diffusiveFlux( entity, quadrature[ pt ], u_x, diffusiveFlux );
        
        // Multiply by all base function gradients
        for( unsigned int i = 0; i < numDofs; ++i )
        {
          // evaluate gradient of base function
          RangeJacobianType phi_x;
          baseFunctionSet.jacobian( i, quadrature[ pt ], phi_x );
          phi_x[ 0 ] = FMatrixHelp :: mult( inv, phi_x[ 0 ] );
         
          // update destination function
          w[ i ] += weight * (diffusiveFlux[ 0 ] * phi_x[ 0 ]);
        }
      }
    }

    template< class EntityType, class LocalMatrixType >
    inline void assembleMatrix ( const EntityType &entity,
                                 LocalMatrixType &localMatrix ) const
    {
      // type of base function sets
      typedef typename LocalMatrixType :: DomainBaseFunctionSetType
        DomainBaseFunctionSetType;
      typedef typename LocalMatrixType :: RangeBaseFunctionSetType
        RangeBaseFunctionSetType;
      typedef typename LocalMatrixType :: RangeFieldType FieldType;

      // type of jacobians
      typedef typename DomainBaseFunctionSetType :: JacobianRangeType DomainJacobianType;
      typedef typename RangeBaseFunctionSetType :: JacobianRangeType RangeJacobianType;

      // obtain geometry from the entity
      const GeometryType &geometry = entity.geometry();

      // obtain base function set
      const DomainBaseFunctionSetType &domainBaseFunctionSet
        = localMatrix.domainBaseFunctionSet();
      const RangeBaseFunctionSetType &rangeBaseFunctionSet
        = localMatrix.rangeBaseFunctionSet();

      const unsigned int rows = localMatrix.rows();
      const unsigned int columns = localMatrix.columns();

      // clear the locl matrix
      localMatrix.clear();

      // Loop filling the local matrix
      QuadratureType quadrature( entity, 2 * polynomialOrder );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // get quadrature point
        const QuadraturePointType &point = quadrature.point( pt );

        // get jacobian inverse of reference mapping
        const GeometryJacobianInverseType &inv
          = geometry.jacobianInverseTransposed( point );

        // weight of this point in the integral
        const FieldType weight
          = geometry.integrationElement( point ) * quadrature.weight( pt );
       
        // Loop over all base functions for matrix rows
        for( unsigned int i = 0; i < rows; ++i )
        {
          // evaluate gradient of i-th base function
          RangeJacobianType phi_x;
          rangeBaseFunctionSet.jacobian( i, quadrature[ pt ], phi_x );
          phi_x[ 0 ] = FMatrixHelp :: mult( inv, phi_x[ 0 ] );
         
          // Loop over all base functions for matrix columns
          for( unsigned int j = 0; j < columns; ++j )
          {
            // evaluate gradient of j-th base function
            DomainJacobianType psi_x;
            domainBaseFunctionSet.jacobian( j, quadrature[ pt ], psi_x );
            
            // calculate diffusive flux
            RangeJacobianType diffusiveFlux;
            model_.diffusiveFlux( entity, quadrature[ pt ], psi_x, diffusiveFlux );
            diffusiveFlux[ 0 ] = FMatrixHelp :: mult( inv, diffusiveFlux[ 0 ] );
 
            // update the local matrix
            localMatrix.add( i, j, weight * (diffusiveFlux[ 0 ] * phi_x[ 0 ]) );
          }
        }
      }
    }
  };



  template< class TraitsImp, class ModelImp >
  class DiffusionOperator
  : public IntegrationOperator
    < DefaultIntegrationOperatorTraits
      < LocalDiffusionOperator< TraitsImp, ModelImp >,
        typename TraitsImp :: DiscreteFunctionType
      >,
      true
    >
  {
  public:
    typedef TraitsImp TraitsType;
    typedef ModelImp ModelType;

    typedef LocalDiffusionOperator< TraitsType, ModelType> LocalOperatorType;

    typedef typename TraitsType :: DiscreteFunctionType DiscreteFunctionType;

  private:
    typedef DefaultIntegrationOperatorTraits< LocalOperatorType, DiscreteFunctionType >
      IntegrationOperatorTraitsType;
 
    typedef DiffusionOperator< TraitsType, ModelType > ThisType;
    typedef IntegrationOperator< IntegrationOperatorTraitsType, true > BaseType;

  public:
    typedef typename TraitsType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType DomainFunctionSpaceType;
    typedef DiscreteFunctionSpaceType RangeFunctionSpaceType;

    typedef DiscreteFunctionType DomainFunctionType;
    typedef DiscreteFunctionType RangeFunctionType;

    typedef typename TraitsType :: SourceFunctionType SourceFunctionType;

    typedef IdentityOperator< DiscreteFunctionType > DomainProjectionType;
    typedef EllipticSourceProjection
      < DefaultEllipticSourceProjectionTraits< SourceFunctionType, DiscreteFunctionType > >
      RangeProjectionType;

  protected:
    LocalOperatorType localOperator_;

  public:
    inline DiffusionOperator ( const DiscreteFunctionSpaceType &dfSpace,
                               const ModelType &model )
    : BaseType( localOperator_, dfSpace, dfSpace ),
      localOperator_( model )
    {
    }

    inline const DomainProjectionType domainProjection () const
    {
      return DomainProjectionType();
    }

    inline const RangeProjectionType rangeProjection () const
    {
      return RangeProjectionType();
    }
  };
  
}

#endif
