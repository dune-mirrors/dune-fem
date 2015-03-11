#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#error "This file seems to be unused and will be removed!!!"

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

// EllipticOperator
// ----------------

template< class DiscreteFunction, class Model >
struct EllipticOperator
: public Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType::Entity::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry IntersectionGeometryType;

  typedef Dune::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::CachingQuadrature< GridPartType, 0 > QuadratureType;

public:
  explicit EllipticOperator ( const ModelType &model = Model() )
  : model_( model )
  {}

  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  template< class JacobianOperator >
  void jacobian ( const DiscreteFunction &u, JacobianOperator &jOp ) const;
private:
  ModelType model_;
  double beta_;
};

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
{
  w.clear();

  const DiscreteFunctionSpaceType &dfSpace = w.space();
  for( const auto& entity : dfSpace )
  {
    const GeometryType &geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalFunctionType wLocal = w.localFunction( entity );

    const int quadOrder = uLocal.order() + wLocal.order();

    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType vu, avu;
        uLocal.evaluate( quadrature[ pt ], vu );
        model_.massFlux( entity, quadrature[ pt ], vu, avu );

        JacobianRangeType du, adu;
        uLocal.jacobian( quadrature[ pt ], du );
        model_.diffusiveFlux( entity, quadrature[ pt ], du, adu );

        avu *= weight;
        adu *= weight;
        wLocal.axpy( quadrature[ pt ], avu, adu );
      }
    }
  }
  w.communicate();
}

// MassModel
// --------------

template< class FunctionSpace >
struct MassModel
{
  typedef FunctionSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  template< class Entity, class Point >
  void massFlux ( const Entity &entity, const Point &x,
                  const RangeType &value, RangeType &flux ) const
  {
    flux = value;
  }

  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity, const Point &x,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    flux = 0;
  }
};


template< class Function, class DiscreteFunction >
void assembleRHS ( const Function &function, DiscreteFunction &rhs )
{
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::IteratorType::Entity::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::CachingQuadrature< GridPartType, 0 > QuadratureType;

  rhs.clear();

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();
  for( const auto& entity : dfSpace )
  {
    const GeometryType &geometry = entity.geometry();

    typename Function::LocalFunctionType localFunction =
             function.localFunction( entity);
    LocalFunctionType rhsLocal = rhs.localFunction( entity );

    QuadratureType quadrature( entity, 2*dfSpace.order()+1 );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      // obtain quadrature point
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

      // evaluate f
      typename Function::RangeType f;
      localFunction.evaluate( quadrature[ pt ], f );

      // multiply by quadrature weight
      f *= quadrature.weight( pt ) * geometry.integrationElement( x );

      // add f * phi_i to rhsLocal[ i ]
      rhsLocal.axpy( quadrature[ pt ], f );
    }
  }
  rhs.communicate();
}
template< class DiscreteFunction, class Model >
template< class JacobianOperator >
void EllipticOperator< DiscreteFunction, Model >
  ::jacobian ( const DiscreteFunction &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  const DiscreteFunctionSpaceType &dfSpace = u.space();

  jOp.reserve();
  jOp.clear();

  std::vector< typename LocalFunctionType::RangeType > phi( dfSpace.mapper().maxNumDofs() );
  std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace.mapper().maxNumDofs() );

  for( const auto& entity : dfSpace )
  {
    const GeometryType &geometry = entity.geometry();

    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const BaseFunctionSetType &baseSet = jLocal.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      const typename GeometryType::Jacobian &gjit = geometry.jacobianInverseTransposed( x );

      // evaluate all basis functions at given quadrature point
      baseSet.evaluateAll( quadrature[ pt ], phi );

      // evaluate jacobians of all basis functions at given quadrature point
      baseSet.jacobianAll( quadrature[ pt ], gjit, dphi );

      for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
      {
        typename LocalFunctionType::RangeType aphi;
        model_.massFlux( entity, quadrature[ pt ], phi[ localCol ], aphi );
        typename LocalFunctionType::JacobianRangeType adphi;
        model_.diffusiveFlux( entity, quadrature[ pt ], dphi[ localCol ], adphi );

        // get column object and call axpy method
        jLocal.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
      }
    }
  }
}

#endif // #ifndef ELLIPTIC_HH
