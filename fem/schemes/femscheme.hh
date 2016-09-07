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

#include <iostream>
#include <memory>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/schemes/solver.hh>

// estimator for residual
#include "estimator.hh"

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem/schemes/dgrhs.hh>
#include <dune/fem/schemes/diffusionmodel.hh>

// FemScheme
//----------

/*******************************************************************************
 * template arguments are:
 * - GridPart: the part of the grid used to tesselate the
 *             computational domain
 * - Model: description of the data functions and methods required for the
 *          elliptic operator (massFlux, diffusionFlux)
 *******************************************************************************/
template< class Space, class Model,
  template<class LinOp,class M,class Constraint = Dune::DirichletConstraints< M, typename LinOp::RangeFunctionType::DiscreteFunctionSpaceType> > class DifferentiableOperator,
    SolverType solver >
class FemScheme
{
public:
  //! type of the mathematical model
  typedef Model ModelType;
  typedef typename ModelType::ExactSolutionType ExactSolutionType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;
  static_assert( std::is_same< typename Space::GridPartType, GridPartType >::value,
        "GridPart of Space has to be identical to GridPart of Model class" );

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename ModelType::FunctionSpaceType FunctionSpaceType;

  //! choose type of discrete function space
  typedef Space DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Solvers<DiscreteFunctionSpaceType,solver,false> UsedSolverType;
  static_assert( UsedSolverType::solverConfigured, "chosen solver is not configured" );

  typedef typename UsedSolverType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename UsedSolverType::LinearOperatorType LinearOperatorType;

  //! type of restriction/prolongation projection for adaptive simulations
  //! (use default here, i.e. LagrangeInterpolation)
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  //! type of error estimator
  typedef Estimator< DiscreteFunctionType, Model, GridPartType::dimension == GridPartType::dimensionworld > EstimatorType;

  static const int dimRange = FunctionSpaceType::dimRange;
  typedef DiscreteFunctionType SolutionType;
  /*********************************************************/

  FemScheme ( const DiscreteFunctionSpaceType &space, const ModelType &model, const std::string &name, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : model_( model ),
    name_( name ),
    space_( space ),
    rhs_( "rhs", space_ ),
    // the elliptic operator (implicit)
    implicitOperator_( new DifferentiableOperator< LinearOperatorType, ModelType >( model_, space_ ) ),
    // create linear operator (domainSpace,rangeSpace)
    linearOperator_( new LinearOperatorType( "assembled elliptic operator", space_, space_) ), // , parameter ) ),
    estimator_( space_, model ),
    exactSolution_( model_.exactSolution( gridPart() ) ),
    parameter_(parameter)
  {}

  const ExactSolutionType &exactSolution() const { return exactSolution_; }

  void prepare ()
  {
    typedef DifferentiableOperator< LinearOperatorType, ModelType > OperatorType;
    // assemble rhs
    assembleDGRHS( model_, model_.rightHandSide( gridPart() ), model_.neumanBoundary( gridPart() ), rhs_, 20 );

    // set boundary values to the rhs - since implicitOperator is of
    // abstract base type we need to cast here
    dynamic_cast< OperatorType & >( *implicitOperator_ ).prepare( rhs_ );
  }

  void prepare ( const DiscreteFunctionType &add )
  {
    prepare();
    rhs_ += add;
  }

  void operator() ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) { (*implicitOperator_)( arg, dest ); }

  void solve ( DiscreteFunctionType &solution, bool assemble )
  {
    typedef DifferentiableOperator< LinearOperatorType, ModelType > OperatorType;
    typedef typename UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
    typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > InverseOperatorType;
    InverseOperatorType invOp( dynamic_cast< OperatorType & >( *implicitOperator_ ), parameter_ );
    invOp( rhs_, solution );
  }

  //! mark elements for adaptation
  bool mark ( double tolerance ) { return estimator_.mark( tolerance ); }

  //! calculate error estimator
  double estimate ( const DiscreteFunctionType &solution ) { return estimator_.estimate( solution ); }

  const std::string &name () { return name_; }

  const GridPartType &gridPart () const { return space().gridPart(); }
  const DiscreteFunctionSpaceType &space( ) const { return space_; }

protected:
  const ModelType &model_;   // the mathematical model
  const std::string name_;
  const DiscreteFunctionSpaceType &space_; // discrete function space
  DiscreteFunctionType rhs_;        // the right hand side
  std::unique_ptr< Dune::Fem::DifferentiableOperator< LinearOperatorType > > implicitOperator_;
  std::unique_ptr< Dune::Fem::Operator< DiscreteFunctionType,DiscreteFunctionType > > linearOperator_;
  EstimatorType estimator_; // estimator for residual error
  const ExactSolutionType exactSolution_;
  const Dune::Fem::ParameterReader parameter_;
};

#endif // #ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
