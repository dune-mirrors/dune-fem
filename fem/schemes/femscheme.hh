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

// estimator for residual <- not adapted to new diffusion model yet
#include "estimator.hh"

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem/schemes/diffusionmodel.hh>

// FemScheme
//----------

template< class Operator, SolverType s >
class FemScheme
{
public:
  static const SolverType solver = s;
  //! type of the mathematical model
  typedef typename Operator::ModelType ModelType;
  typedef typename Operator::RangeDiscreteFunctionType DiscreteFunctionType;
  typedef Operator DifferentiableOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;
  static_assert( std::is_same< typename DiscreteFunctionSpaceType::GridPartType, GridPartType >::value,
        "GridPart of Space has to be identical to GridPart of Model class" );

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Solvers<DiscreteFunctionSpaceType,solver,false> UsedSolverType;
  static_assert( UsedSolverType::solverConfigured, "chosen solver is not configured" );

  typedef typename UsedSolverType::LinearOperatorType LinearOperatorType;

  //! type of restriction/prolongation projection for adaptive simulations
  //! (use default here, i.e. LagrangeInterpolation)
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  //! type of error estimator
  typedef Estimator< DiscreteFunctionType, ModelType, GridPartType::dimension == GridPartType::dimensionworld > EstimatorType;

  static const int dimRange = FunctionSpaceType::dimRange;
  typedef DiscreteFunctionType SolutionType;
  /*********************************************************/

  FemScheme ( const DiscreteFunctionSpaceType &space, const ModelType &model, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : model_( model ),
    space_( space ),
    // the elliptic operator (implicit)
    implicitOperator_( model_, space_ ),
    // create linear operator (domainSpace,rangeSpace)
    linearOperator_( "assembled elliptic operator", space_, space_ ), // , parameter ),
    estimator_( space_, model ),
    parameter_(parameter)
  {}

  const DifferentiableOperatorType &fullOperator() const
  {
    return implicitOperator_;
  }

  void constraint( DiscreteFunctionType &u ) const
  {
    implicitOperator_.prepare( u );
  }

  void operator() ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
  {
    implicitOperator_( arg, dest );
  }
  template <class GridFunction>
  void operator() ( const GridFunction &arg, DiscreteFunctionType &dest ) const
  {
    implicitOperator_.apply( arg, dest );
  }

  void solve ( DiscreteFunctionType &solution ) const
  {
    typedef typename UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
    typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > InverseOperatorType;
    InverseOperatorType invOp( implicitOperator_, parameter_ );
    DiscreteFunctionType bnd(solution);
    bnd.clear();
    implicitOperator_.prepare( bnd );
    invOp( bnd, solution );
  }

  template <class GridFunction>
  const LinearOperatorType &assemble( const GridFunction &ubar )
  {
    implicitOperator_.apply(ubar, linearOperator_);
    return linearOperator_;
  }

  //! mark elements for adaptation
  bool mark ( double tolerance ) { return estimator_.mark( tolerance ); }

  //! calculate error estimator
  double estimate ( const DiscreteFunctionType &solution ) { return estimator_.estimate( solution ); }

  const GridPartType &gridPart () const { return space().gridPart(); }
  const DiscreteFunctionSpaceType &space( ) const { return space_; }

protected:
  const ModelType &model_;   // the mathematical model
  const DiscreteFunctionSpaceType &space_; // discrete function space
  DifferentiableOperatorType implicitOperator_;
  LinearOperatorType linearOperator_;
  EstimatorType estimator_; // estimator for residual error
  const Dune::Fem::ParameterReader parameter_;
};

#endif // #ifndef DUNE_FEM_SCHEMES_FEMSCHEME_HH
