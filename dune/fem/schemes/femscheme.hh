/**************************************************************************

  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2015 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger
  Copyright (C) 2004 - 2015 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2015 Mirko Kraenkel
  Copyright (C) 2006 - 2015 Christoph Gersbacher
  Copyright (C) 2006 - 2015 Martin Nolte
  Copyright (C) 2011 - 2015 Tobias Malkmus
  Copyright (C) 2012 - 2015 Stefan Girke
  Copyright (C) 2013 - 2015 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner
  Copyright (C) 2015        Marco Agnese
  Copyright (C) 2015        Martin Alkaemper


  The dune-fem module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

**************************************************************************/
#ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
#define DUNE_FEM_SCHEMES_FEMSCHEME_HH

#include <type_traits>
#include <utility>
#include <iostream>
#include <memory>
#include <dune/common/typeutilities.hh>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

// FemScheme
//----------

template < class Op, class DF, typename = void >
struct AddDirichletBC
{
  static const bool value = false;
  using DirichletBlockVector = void;
};
template < class Op, class DF>
struct AddDirichletBC<Op,DF,std::enable_if_t<std::is_void< decltype( std::declval<const Op>().
            setConstraints( std::declval<DF&>() ) )>::value > >
{
  static const bool value = true;
  using DirichletBlockVector = typename Op::DirichletBlockVector;
};


template< class Operator, class LinearInverseOperator >
class FemScheme
{
public:
  //! type of the mathematical model
  typedef typename Operator::ModelType ModelType;
  typedef typename Operator::DomainFunctionType DomainFunctionType;
  typedef typename Operator::RangeFunctionType  RangeFunctionType;
  typedef typename Operator::RangeFunctionType  DiscreteFunctionType;
  typedef Operator DifferentiableOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef LinearInverseOperator LinearInverseOperatorType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;
  static_assert( std::is_same< typename DiscreteFunctionSpaceType::GridPartType, GridPartType >::value,
        "GridPart of Space has to be identical to GridPart of Model class" );

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename Operator::JacobianOperatorType JacobianOperatorType;
  typedef typename Operator::JacobianOperatorType LinearOperatorType;
  typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > InverseOperatorType;
  typedef typename InverseOperatorType::ErrorMeasureType ErrorMeasureType;

  typedef typename FunctionSpaceType::RangeType RangeType;
  static const int dimRange = FunctionSpaceType::dimRange;
  static constexpr bool addDirichletBC = AddDirichletBC<Operator,DomainFunctionType>::value;
  using DirichletBlockVector = typename AddDirichletBC<Operator,DomainFunctionType>::DirichletBlockVector;
  /*********************************************************/

  FemScheme ( const DiscreteFunctionSpaceType &space, ModelType &model, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : space_( space ),
    // the elliptic operator (implicit)
    implicitOperator_( space, space, model, parameter ),
    // create linear operator (domainSpace,rangeSpace)
    invOp_( parameter )
  {}

  const DifferentiableOperatorType &fullOperator() const { return implicitOperator_; }
  DifferentiableOperatorType &fullOperator() { return implicitOperator_; }

  template <typename O = DifferentiableOperatorType>
  auto setQuadratureOrders(unsigned int interior, unsigned int surface)
  -> Dune::void_t< decltype( std::declval< O >().setQuadratureOrders(0,0) ) >
  {
    fullOperator().setQuadratureOrders(interior,surface);
  }

  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  setConstraints( DomainFunctionType &u ) const
  {
    implicitOperator_.setConstraints( u );
  }
  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
  {
    implicitOperator_.setConstraints( u, v );
  }
  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
  {
    implicitOperator_.subConstraints( u, v );
  }
  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  setConstraints( const RangeType &value, DiscreteFunctionType &u ) const
  {
    implicitOperator_.setConstraints( value, u );
  }
  // template <typename O = Operator>
  // std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,
  //    const decltype(implicitOperator_.dirichletBlocks())&>
  const auto &dirichletBlocks() const
  {
    if constexpr (AddDirichletBC<Operator,DomainFunctionType>::value)
      return implicitOperator_.dirichletBlocks();
  }

  void operator() ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
  {
    implicitOperator_( arg, dest );
  }
  template <class GridFunction>
  auto operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
  -> Dune::void_t<decltype(std::declval<const Operator&>()(arg,dest))>
  {
    implicitOperator_( arg, dest );
  }

  struct SolverInfo
  {
    SolverInfo(bool pconverged,int plinearIterations,int pnonlinearIterations)
      : converged(pconverged), linearIterations(plinearIterations), nonlinearIterations(pnonlinearIterations)
    {}
    bool converged;
    int linearIterations;
    int nonlinearIterations;
  };
  void setErrorMeasure(ErrorMeasureType &errorMeasure) const
  {
    invOp_.setErrorMeasure(errorMeasure);
  }
  SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution ) const
  {
    invOp_.bind( implicitOperator_ );
    DiscreteFunctionType rhs0 = rhs;
    setZeroConstraints( rhs0 );
    setModelConstraints( solution );
    invOp_( rhs0, solution );
    invOp_.unbind();
    return SolverInfo(invOp_.converged(),invOp_.linearIterations(),invOp_.iterations());
  }
  SolverInfo solve ( DiscreteFunctionType &solution ) const
  {
    DiscreteFunctionType bnd( solution );
    bnd.clear();
    setModelConstraints( solution );
    invOp_.bind(fullOperator());
    invOp_( bnd, solution );
    invOp_.unbind();
    return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations() );
  }

  template< class GridFunction, std::enable_if_t<
        std::is_same< decltype(
          std::declval< const DifferentiableOperatorType >().jacobian(
              std::declval< const GridFunction& >(), std::declval< JacobianOperatorType& >()
            )
          ), void >::value, int> i = 0
    >
  void jacobian( const GridFunction &ubar, JacobianOperatorType &linOp ) const
  {
    implicitOperator_.jacobian(ubar, linOp);
  }

  const GridPartType &gridPart () const { return space().gridPart(); }
  const DiscreteFunctionSpaceType &space( ) const { return space_; }

  const ModelType &model() const
  {
    return implicitOperator_.model();
  }
  ModelType &model()
  {
    return implicitOperator_.model();
  }
protected:
  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  setZeroConstraints( DiscreteFunctionType &u ) const { implicitOperator_.setConstraints( RangeType(0), u ); }
  template<class...Args>
  void setZeroConstraints(Args&&...) const { }
  template <typename O = Operator>
  std::enable_if_t<AddDirichletBC<O,DomainFunctionType>::value,void>
  setModelConstraints( DiscreteFunctionType &u ) const { fullOperator().setConstraints( u ); }
  template<class...Args>
  void setModelConstraints(Args&&... ) const { }
  const DiscreteFunctionSpaceType &space_; // discrete function space
  DifferentiableOperatorType implicitOperator_;
  mutable InverseOperatorType invOp_;
};

#endif // #ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
