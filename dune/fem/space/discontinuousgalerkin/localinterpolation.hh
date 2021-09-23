#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH

// dune-fem includes
#include <dune/grid/common/capabilities.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/agglomerationquadrature.hh>

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
    template< class DiscreteFunctionSpace >
    class DiscontinuousGalerkinLocalInterpolation
    {
      typedef DiscontinuousGalerkinLocalInterpolation< DiscreteFunctionSpace > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridType  GridType;
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

      static const bool isAlwaysAffine = Dune::Capabilities::isCartesian< GridType >::v ||
         ( Dune::Capabilities::hasSingleGeometryType< GridType >::v &&  ((Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId >> 1) == 0)) ;
      // always true for orthonormal spaces
      //static const bool isAlwaysAffine = true;

    private:
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename RangeType::value_type RangeFieldType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType :: LocalMassMatrixType
        LocalMassMatrixType;

      typedef typename LocalMassMatrixType::VolumeQuadratureType  QuadratureType;

    public:
      DiscontinuousGalerkinLocalInterpolation ( const DiscreteFunctionSpaceType &space )
      : space_( space )
      {}

      DiscontinuousGalerkinLocalInterpolation ( const ThisType &other ) = default;
      DiscontinuousGalerkinLocalInterpolation ( ThisType &&other ) = default;

      void bind( const EntityType& ) {}
      void unbind() {}

      ThisType &operator= ( const ThisType &other ) = delete;

      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &localFunction, LocalDofVector &dofs ) const
      {
        // set all dofs to zero
        std::fill( dofs.begin(), dofs.end(), typename LocalDofVector::value_type(0) );

        // get entity and geometry
        const EntityType &entity = localFunction.entity();

        if( entity.type().isNone() )
        {
          typedef ElementQuadrature< GridPartType, EntityType::codimension > ElementQuadratureType;
          ElementQuadratureType quadrature( entity, localFunction.order() + space_.order( entity ) );
          bool isAffine = computeInterpolation( entity, quadrature, localFunction, dofs );
          if( ! isAffine )
          {
            typedef LocalMassMatrix< DiscreteFunctionSpaceType, ElementQuadratureType > AggloMassMatrix;
            AggloMassMatrix massMat( space_, massMatrix().volumeQuadratureOrder( entity ) );
            // apply inverse of mass matrix
            auto basisFunctionSet = space_.basisFunctionSet(entity);
            massMat.applyInverse( entity, basisFunctionSet, dofs );
          }
        }
        else
        {
          QuadratureType quadrature( entity, localFunction.order() + space_.order( entity ) );
          bool isAffine = computeInterpolation( entity, quadrature, localFunction, dofs );
          if( ! isAffine )
          {
            // apply inverse of mass matrix
            auto basisFunctionSet = space_.basisFunctionSet(entity);
            massMatrix().applyInverse( entity, basisFunctionSet, dofs );
          }
        }

      }

    private:
      template<class QuadImpl, class LocalFunction, class LocalDofVector >
      bool computeInterpolation( const EntityType& entity,
                                 const QuadImpl& quadrature,
                                 const LocalFunction &localFunction,
                                 LocalDofVector &dofs ) const
      {
        const int nop = quadrature.nop();

        auto& values = space_.localMassMatrixStorage().second;

        // adjust size of values
        values.resize( nop );

        // evaluate local function for all quadrature points
        localFunction.evaluateQuadrature( quadrature, values );

        bool isAffine = isAlwaysAffine ;
        if( ! isAlwaysAffine )
        {
          const auto geometry = entity.geometry();
          isAffine = geometry.affine();

          if( ! isAffine )
          {
            // apply weight
            for(auto qp : quadrature )
              values[ qp.index() ] *= qp.weight() * geometry.integrationElement( qp.position() );
          }
        }

        if( isAffine )
        {
          // apply weight only
          for(auto qp : quadrature )
            values[ qp.index() ] *= qp.weight();
        }

        // add values to local function
        space_.basisFunctionSet(entity).axpy( quadrature, values, dofs );

        return isAffine;
      }

      const LocalMassMatrixType &massMatrix () const
      {
        return space_.localMassMatrixStorage().first;
      }

      const DiscreteFunctionSpaceType& space_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALINTERPOLATION_HH
