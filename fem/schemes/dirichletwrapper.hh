#ifndef DUNE_FEM_DIRICHLETWRAPPER_HH
#define DUNE_FEM_DIRICHLETWRAPPER_HH

#include <cstddef>

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/schemes/dirichletconstraints.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>


template< class Operator >
struct DirichletWrapperOperator
: public Dune::Fem::DifferentiableOperator< typename Operator::JacobianOperatorType >
{
  typedef typename Operator::DomainFunctionType DomainFunctionType;
  typedef typename Operator::RangeFunctionType  RangeFunctionType;
  typedef typename Operator::ModelType ModelType;
  typedef typename Operator::DirichletModelType DirichletModelType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;
  typedef typename Operator::JacobianOperatorType JacobianOperatorType;
  typedef typename RangeDiscreteFunctionSpaceType::RangeType DomainRangeType;
  typedef Dune::DirichletConstraints< ModelType, RangeDiscreteFunctionSpaceType > ConstraintsType;

  template <class... Args>
  DirichletWrapperOperator ( Args&&... args )
    : op_( std::forward<Args&>(args)... ) , constraints_( op_.model(), op_.rangeSpace() )
  {}

  void setConstraints( DomainFunctionType &u ) const
  {
    // set boundary values for solution from model
    constraints()( u );
  }
  void setConstraints( const DomainRangeType &value, DomainFunctionType &u ) const
  {
    // set values for solution to a given constant value
    constraints()( value, u );
  }
  template <class GF>
  void setConstraints( const GF &u, RangeFunctionType &w ) const
  {
    // set boundary values for solution from a general grid function
    constraints()( u, w, ConstraintsType::Operation::set );
  }
  template <class GF>
  void subConstraints( const GF &u, RangeFunctionType &w ) const
  {
    // subtract boundary values from solution
    constraints()( u, w, ConstraintsType::Operation::sub );
  }

  //! application operator
  virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
  {
    op_(u,w);
    subConstraints( u, w );
  }
  template <class GF>
  auto operator()( const GF &u, RangeFunctionType &w ) const
  -> Dune::void_t<decltype(std::declval<const Operator&>()(u,w))>
  {
    op_(u,w);
    subConstraints( u, w );
  }

  void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
  {
    op_.jacobian(u,jOp);
    constraints().applyToOperator( jOp );
  }
  template <class GridFunctionType>
  auto jacobian ( const GridFunctionType &u, JacobianOperatorType &jOp ) const
  -> Dune::void_t<decltype(std::declval<const Operator&>().jacobian(u,jOp))>
  {
    op_.jacobian(u,jOp);
    constraints().applyToOperator( jOp );
  }

  const DomainDiscreteFunctionSpaceType& domainSpace() const
  {
    return op_.domainSpace();
  }
  const RangeDiscreteFunctionSpaceType& rangeSpace() const
  {
    return op_.rangeSpace();
  }

  ModelType &model () const { return op_.model(); }
  const ConstraintsType &constraints () const { return constraints_; }

private:
  Operator op_;
  ConstraintsType constraints_;
};
#endif // #ifndef DUNE_FEM_CONSTRAINTSWRAPPER_HH
