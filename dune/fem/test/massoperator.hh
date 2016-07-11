#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

#include <dune/fem/quadrature/lifted.hh>


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

  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > local( dfSpace_ );
  typename Function::LocalFunctionType uLocal( u );

  w.clear();

  // run over entities
  for( const EntityType &entity : dfSpace_ )
  {
    const GeometryType &geometry = entity.geometry();

    local.init( entity );
    local.clear();

    uLocal.init( entity );

    // run over quadrature points
    QuadratureType quadrature( entity, 2*dfSpace_.order()+1 );
    for( const auto qp : quadrature )
    {
      // evaluate u
      const typename QuadratureType::CoordinateType &x = qp.position();

      RangeType uValue;
      uLocal.evaluate( qp, uValue );

      // put all things together and don't forget quadrature weights
      const FieldType weight = qp.weight()*geometry.integrationElement( x );

      // apply weight
      uValue *= weight;

      // add to local function
      local.axpy( qp, uValue );
    }

    w.addLocalDofs( entity, local.localDofVector() );
  }

  w.communicate();
}


template< class DiscreteFunction, class LinearOperator >
void MassOperator< DiscreteFunction, LinearOperator >::assemble ()
{
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType FieldType;

  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > LocalMatrixType;

  BaseType::reserve( Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>( dfSpace_, dfSpace_ ) );
  BaseType::clear();

  LocalMatrixType localMatrix( dfSpace_, dfSpace_ );

  typename Dune::Fem::IdTensor< typename DiscreteFunctionSpaceType::RangeType > mass;

  // run over entities
  for( const EntityType &entity : dfSpace_ )
  {
    const GeometryType &geometry = entity.geometry();

    localMatrix.init( entity, entity );
    localMatrix.clear();

    // run over quadrature points
    for( const auto &qp : liftedQuadrature( QuadratureType( entity, 2*dfSpace_.order() ), geometry ) )
    {
      // get quadrature weight
      const FieldType weight = qp.weight();

      // update system matrix
      localMatrix.axpy( qp,
          std::make_pair(
            // Dune::Fem::scaledTensor< typename DiscreteFunctionSpaceType::RangeType > ( weight, mass ),
            mass,
            Dune::Fem::ZeroTensor< typename DiscreteFunctionSpaceType::RangeType > () ),
          std::make_pair(
            Dune::Fem::ZeroTensor< typename DiscreteFunctionSpaceType::JacobianRangeType > (),
            Dune::Fem::ZeroTensor< typename DiscreteFunctionSpaceType::JacobianRangeType > () ) );
            // Dune::Fem::scaledTensor< typename DiscreteFunctionSpaceType::JacobianRangeType >( weight, diffusion ) ) );
    }

    // add to global matrix
    BaseType::addLocalMatrix( entity, entity, localMatrix );
  }
  BaseType::communicate();
}

#endif // #ifndef MASSOPERATOR_HH
