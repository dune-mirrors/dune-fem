#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH

// dune-fem includes
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Local interpolation for Discontinuous Galerkin spaces
*/

namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinLocalInterpolation
    // ---------------------------------------

    /**
     * Local interpolation for Discontinuous Galerkin spaces.
     */
    template< class DiscreteFunctionSpace, template< class, int > class Quadrature = CachingQuadrature >
    class DiscontinuousGalerkinLocalInterpolation 
    {
      typedef DiscontinuousGalerkinLocalInterpolation< DiscreteFunctionSpace, Quadrature > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

    private:
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename RangeType::value_type RangeFieldType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef Quadrature< GridPartType, EntityType::codimension > QuadratureType;

      typedef LocalMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

    public:
      DiscontinuousGalerkinLocalInterpolation ( const DiscreteFunctionSpaceType &space, const int order = -1 )
      : order_( order < 0 ? 2*space.order() : order ),
        massMatrix_( space, order_ )
      {}

    private:
      // forbid copying
      DiscontinuousGalerkinLocalInterpolation ( const ThisType &other );
      // forbid assignment
      ThisType &operator= ( const ThisType &other );

    public:
      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        // get entity and geometry
        const EntityType &entity = localFunction.entity();
        const typename EntityType::Geometry geometry = entity.geometry();
        
        QuadratureType quadrature( entity, order_ );
        const int nop = quadrature.nop();
        for( int qp = 0; qp < nop; ++qp )
        {
          // evaluate local function
          RangeType value;
          localFunction.evaluate( quadrature[ qp ], value );

          // compute quadrature weight
          const RangeFieldType intel = quadrature.weight( qp )*geometry.integrationElement( quadrature.point( qp ) );
          value *= intel;

          dofs.axpy( quadrature[ qp ], value );
        }

        // apply inverse of mass matrix
        massMatrix().applyInverse( entity, dofs );
      }

    private:
      const LocalMassMatrixType &massMatrix () const { return massMatrix_; }

      const int order_;
      LocalMassMatrixType massMatrix_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH
