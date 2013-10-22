#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#if defined HAVE_PETSC

#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>


template< class DiscreteFunction >
class MassOperator
: public Dune::Fem:: PetscLinearOperator< DiscreteFunction, DiscreteFunction >
{
  typedef MassOperator< DiscreteFunction > ThisType;
  typedef Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction > BaseType;

  typedef DiscreteFunction DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef Dune::Fem::CachingQuadrature< typename DiscreteFunctionSpaceType::GridPartType, 0 > QuadratureType;
public:
  explicit MassOperator ( const DiscreteFunctionSpaceType &dfSpace )
  : BaseType( "mass operator", dfSpace , dfSpace ),
    dfSpace_( dfSpace )
  {
    assemble();
  }

  template< class Function >
  void assembleRHS( const Function &u, DiscreteFunctionType &w ) const;

private:
  void assemble (); 

  const DiscreteFunctionSpaceType &dfSpace_;
};


template< class DiscreteFunction >
template< class Function >
void MassOperator< DiscreteFunction >
  ::assembleRHS( const Function &u, DiscreteFunctionType &w ) const
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  w.clear();

  // run over entities
  const IteratorType end = dfSpace_.end();
  for( IteratorType it = dfSpace_.begin(); it != end; ++it )
  {

    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalFunctionType localFunction = w.localFunction( entity );

    // run over quadrature points
    QuadratureType quadrature( entity, 2*dfSpace_.order()+1 );
    const unsigned int numQuadraturePoints = quadrature.nop();
    for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      // evaluate u
      const typename QuadratureType::CoordinateType &x = quadrature.point( qp );        

      RangeType uValue;
      u.evaluate( geometry.global( x ), uValue );

      // put all things together and don't forget quadrature weights
      const FieldType weight = quadrature.weight( qp )*geometry.integrationElement( x );

      // apply weight 
      uValue *= weight; 

      // add to local function 
      localFunction.axpy( quadrature[ qp ], uValue );
    }

  }

  w.communicate();
}


template< class DiscreteFunction >
void MassOperator< DiscreteFunction >::assemble ()
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;

  typedef typename BaseType::LocalMatrixType LocalMatrixType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  BaseType::reserve( Dune::Fem::SimpleStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>() );
  BaseType::clear();

  std::vector< typename DiscreteFunctionSpaceType::RangeType > values;

  // run over entities
  const IteratorType end = dfSpace_.end();
  for( IteratorType it = dfSpace_.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalMatrixType localMatrix = BaseType::localMatrix( entity, entity );

    const BasisFunctionSetType &basis = localMatrix.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = basis.size();
    values.resize( numBasisFunctions );

    // run over quadrature points
    QuadratureType quadrature( entity, 2*dfSpace_.order() );
    const unsigned int numQuadraturePoints = quadrature.nop();
    for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      // evaluate base functions
      basis.evaluateAll( quadrature[ qp ], values );

      // get quadrature weight
      const typename QuadratureType::CoordinateType &x = quadrature.point( qp );
      const FieldType weight = quadrature.weight( qp )
                                  * geometry.integrationElement( x );

      // update system matrix
      for( unsigned int i = 0; i < numBasisFunctions; ++i )
      {
        RangeType value = values[ i ];
        // add column
        localMatrix.column( i ).axpy( values, value, weight ); 
      }
    }
  }
  BaseType::communicate();
}

#endif // #if defined HAVE_PETSC

#endif // #ifndef MASSOPERATOR_HH