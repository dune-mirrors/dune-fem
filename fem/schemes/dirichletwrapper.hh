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

#if 0
struct NoConstraints
{
  enum Operation { set = 0, sub = 1 };

  template <class ModelType, class DiscreteFunctionSpaceType>
  NoConstraints( const ModelType&, const DiscreteFunctionSpaceType& )
  {}

  template < class DiscreteFunctionType >
  void operator ()( const DiscreteFunctionType& u) const
  {}

  template < class DiscreteFunctionType >
  void operator ()( const typename DiscreteFunctionType::RangeType& value, DiscreteFunctionType& w ) const
  {}

  template < class GridFunctionType, class DiscreteFunctionType >
  void operator ()( const GridFunctionType& u, DiscreteFunctionType& w, Operation op ) const
  {}

  template <class LinearOperator>
  void applyToOperator( LinearOperator& linearOperator ) const
  {}
};

template< class DomainDiscreteFunction, class RangeDiscreteFunction, class Model>
using DirichletBC = Dune::DirichletConstraints< Model, typename DomainDiscreteFunction::DiscreteFunctionSpaceType >;
template< class DomainDiscreteFunction, class RangeDiscreteFunction, class Model>
using NoBC = NoConstraints;
template< class DomainDiscreteFunction, class RangeDiscreteFunction, class Model>
using DefaultBC = std::conditional_t<
    Model::dimD >= DomainDiscreteFunction::DiscreteFunctionSpaceType::localBlockSize
    && std::is_same<DomainDiscreteFunction,RangeDiscreteFunction>::value,
      DirichletBC<DomainDiscreteFunction,RangeDiscreteFunction,Model>,
      NoBC<DomainDiscreteFunction,RangeDiscreteFunction,Model> >;
#endif


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

  // prepare the solution vector
  void setConstraints( DomainFunctionType &u ) const
  {
    // set boundary values for solution
    constraints()( u );
  }
  // prepare the solution vector
  void setConstraints( const DomainRangeType &value, DomainFunctionType &u ) const
  {
    // set boundary values for solution
    constraints()( value, u );
  }
  template <class GF>
  void setConstraints( const GF &u, RangeFunctionType &w ) const
  {
    // set boundary values for solution
    constraints()( u, w, ConstraintsType::Operation::set );
  }
  template <class GF>
  void subConstraints( const GF &u, RangeFunctionType &w ) const
  {
    // set boundary values for solution
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
