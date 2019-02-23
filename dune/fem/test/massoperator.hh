#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

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
    assemble(*this);
  }

  template< class Function >
  void assembleRHS( const Function &u, DiscreteFunctionType &w ) const;
  void assemble (LinearOperator &matrix) const;
private:
  const DiscreteFunctionSpaceType &dfSpace_;
};

template< class DiscreteFunction, class LinearOperator >
class AffineMassOperator
: public Dune::Fem::DifferentiableOperator<LinearOperator>
{
  public:
  typedef DiscreteFunction DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef LinearOperator JacobianOperatorType;
  template< class Function >
  explicit AffineMassOperator ( const DiscreteFunctionSpaceType &dfSpace, const Function &u )
    : matrix_(dfSpace), rhs_("rhs", dfSpace)
  {
    matrix_.assembleRHS(u, rhs_);
  }
  const DiscreteFunctionType rhs() const { return rhs_; }
  virtual void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
  {
    matrix_(u,w);
    w -= rhs_;
  }
  virtual void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
  {
    matrix_.assemble( jOp );
  }
  private:
  MassOperator<DiscreteFunction,LinearOperator> matrix_;
  DiscreteFunction rhs_;
};

template< class DiscreteFunction, class LinearOperator >
template< class Function >
void MassOperator< DiscreteFunction, LinearOperator >
  ::assembleRHS( const Function &u, DiscreteFunctionType &w ) const
{
  // clear result
  w.clear();

  Dune::Fem::ConstLocalFunction< Function > uLocal( u );
  Dune::Fem::AddLocalContribution< DiscreteFunctionType > wLocal( w );

  // run over entities
  for( const auto &entity : elements( dfSpace_.gridPart(), Dune::Partitions::interiorBorder ) )
  {
    const auto geometry = entity.geometry();

    auto guard = bindGuard( std::tie( uLocal, wLocal ), entity );

    // run over quadrature points
    for( const auto qp : QuadratureType( entity, uLocal.order()+wLocal.order() ) )
    {
      // evaluate u
      RangeType uValue = uLocal.evaluate( qp );

      // apply quadrature weight
      uValue *= qp.weight() * geometry.integrationElement( qp.position() );

      // add to local function
      wLocal.axpy( qp, uValue );
    }
  }
}


template< class DiscreteFunction, class LinearOperator >
void MassOperator< DiscreteFunction, LinearOperator >::assemble (LinearOperator &matrix) const
{
  matrix.reserve( Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>( dfSpace_, dfSpace_ ) );
  matrix.clear();

  Dune::Fem::AddLocalContribution< LinearOperator > localMatrix( matrix );
  Dune::DynamicVector< typename DiscreteFunctionSpaceType::RangeType > values;

  // run over entities
  for( const auto &entity : elements( dfSpace_.gridPart(), Dune::Partitions::interiorBorder ) )
  {
    const auto geometry = entity.geometry();

    auto guard = bindGuard( localMatrix, entity, entity );

    const auto &basis = localMatrix.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = basis.size();
    values.resize( numBasisFunctions );

    // run over quadrature points
    for( const auto qp : QuadratureType( entity, localMatrix.order() ) )
    {
      // evaluate base functions
      basis.evaluateAll( qp, values );

      // apply quadrature weight
      const auto weight = qp.weight() * geometry.integrationElement( qp.position() );
      std::for_each( values.begin(), values.end(), [ weight ] ( auto &a ) { a *= weight; } );

      // update system matrix
      localMatrix.axpy( qp, values );
    }
  }
}

#endif // #ifndef MASSOPERATOR_HH
