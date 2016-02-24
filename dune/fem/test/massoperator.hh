#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>


template< class DiscreteFunction, class LinearOperator >
class MassOperator
: public LinearOperator
{
  typedef MassOperator< DiscreteFunction, LinearOperator > ThisType;
  typedef LinearOperator BaseType;

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


template< class DiscreteFunction, class LinearOperator >
template< class Function >
void MassOperator< DiscreteFunction, LinearOperator >
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
    for( const auto qp : quadrature )
    {
      // evaluate u
      const typename QuadratureType::CoordinateType &x = qp.position();

      RangeType uValue;
      u.evaluate( geometry.global( x ), uValue );

      // put all things together and don't forget quadrature weights
      const FieldType weight = qp.weight()*geometry.integrationElement( x );

      // apply weight
      uValue *= weight;

      // add to local function
      localFunction.axpy( qp, uValue );
    }

  }

  w.communicate();
}


template< class DiscreteFunction, class LinearOperator >
void MassOperator< DiscreteFunction, LinearOperator >::assemble ()
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;

  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > LocalMatrixType;

  BaseType::reserve( Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>( dfSpace_, dfSpace_ ) );
  BaseType::clear();

  LocalMatrixType localMatrix( dfSpace_, dfSpace_ );

  std::vector< typename DiscreteFunctionSpaceType::RangeType > values;

  // run over entities
  const IteratorType end = dfSpace_.end();
  for( IteratorType it = dfSpace_.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    localMatrix.init( entity, entity );
    localMatrix.clear();

    const BasisFunctionSetType &basis = localMatrix.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = basis.size();
    values.resize( numBasisFunctions );

    // run over quadrature points
    QuadratureType quadrature( entity, 2*dfSpace_.order() );
    for( const auto qp : quadrature )
    {
      // evaluate base functions
      basis.evaluateAll( qp, values );

      // get quadrature weight
      const typename QuadratureType::CoordinateType &x = qp.position();
      const FieldType weight = qp.weight()
                                  * geometry.integrationElement( x );

      // update system matrix
      for( unsigned int i = 0; i < numBasisFunctions; ++i )
      {
        RangeType value = values[ i ];
        // add column
        localMatrix.column( i ).axpy( values, value, weight );
      }
    }

    // add to global matrix
    BaseType::addLocalMatrix( entity, entity, localMatrix );
  }
  BaseType::communicate();
}

#endif // #ifndef MASSOPERATOR_HH
