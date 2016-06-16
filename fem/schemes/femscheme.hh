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
#ifndef ELLIPT_FEMSCHEME_HH
#define ELLIPT_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/schemes/solver.hh>

// estimator for residual
#include "estimator.hh"

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem/schemes/elliptic.hh>
#include <dune/fem/schemes/rhs.hh>
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
template < class Space, class Model, int polOrder, SolverType solver >
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

  //! type of adaptation manager handling adapation and DoF compression
  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

  //! type of error estimator
  typedef Estimator< DiscreteFunctionType, Model > EstimatorType;

  static const int dimRange = FunctionSpaceType::dimRange;
  typedef DiscreteFunctionType SolutionType;
  /*********************************************************/

  FemScheme( GridPartType &gridPart,
             const ModelType& implicitModel,
             const std::string &prefix);
  ~FemScheme();
  SolutionType &solution() { return solution_; }
  const DiscreteFunctionType &solution() const { return solution_; }
  const ExactSolutionType& exactSolution() const { return exactSolution_; }

  void prepare();
  void solve ( bool assemble );

  //! mark elements for adaptation
  bool mark ( const double tolerance )
  {
    return estimator_.mark( tolerance );
  }

  //! calculate error estimator
  double estimate()
  {
    return estimator_.estimate( this->implicitModel_.rightHandSide(gridPart_) );
  }

protected:
  const ModelType& implicitModel_;   // the mathematical model
  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with
  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
  DiscreteFunctionType rhs_;        // the right hand side
  Dune::Fem::DifferentiableOperator< LinearOperatorType > *implicitOperator_;
  Dune::Fem::Operator<DiscreteFunctionType,DiscreteFunctionType> *linearOperator_;  // the linear operator (i.e. jacobian of the implicit)
  const double solverEps_ ; // eps for linear solver
  EstimatorType estimator_; // estimator for residual error
  RestrictionProlongationType restrictProlong_ ; // local restriction/prolongation object
  const ExactSolutionType exactSolution_;
};


// DataOutputParameters
// --------------------

template < class Space, class Model, int polOrder, SolverType solver >
FemScheme<Space, Model, polOrder, solver>::
FemScheme( GridPartType &gridPart,
           const ModelType& implicitModel,
           const std::string &prefix)
    : implicitModel_( implicitModel ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( prefix.c_str(), discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_(
          new DifferentiableEllipticOperator< LinearOperatorType, ModelType >
              ( implicitModel_, discreteSpace_ ) ),
      // create linear operator (domainSpace,rangeSpace)
      linearOperator_( new LinearOperatorType( "assembled elliptic operator", discreteSpace_, discreteSpace_ ) ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ),
      estimator_( solution_, implicitModel ),
      restrictProlong_( solution_ ),
      exactSolution_( implicitModel_.exactSolution(gridPart_) )
{
  // set all DoF to zero
  solution_.clear();
}
template < class Space, class Model, int polOrder, SolverType solver >
FemScheme<Space, Model, polOrder, solver>::
~FemScheme()
{
  std::cout << "in FemScheme destructor" << std::endl;
  delete implicitOperator_;
  delete linearOperator_;
}

  //! setup the right hand side
template < class Space, class Model, int polOrder, SolverType solver >
void FemScheme<Space, Model, polOrder, solver>::
prepare()
{
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > OperatorType;
  // assemble rhs
  assembleRHS ( implicitModel_, implicitModel_.rightHandSide(gridPart_), implicitModel_.neumanBoundary(gridPart_), rhs_ );

  // set boundary values to the rhs - since implicitOperator is of
  // abstract base type we need to cast here
  dynamic_cast<OperatorType&>(*implicitOperator_).prepare( implicitModel_.dirichletBoundary(gridPart_), rhs_ );
}

template < class Space, class Model, int polOrder, SolverType solver >
void FemScheme<Space, Model, polOrder, solver>::
solve ( bool assemble )
{
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > OperatorType;
  typedef typename UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
  typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > InverseOperatorType;
  InverseOperatorType invOp( dynamic_cast<OperatorType&>(*implicitOperator_) );
  invOp( rhs_, solution_ );
}

#endif // end #if ELLIPT_FEMSCHEME_HH
