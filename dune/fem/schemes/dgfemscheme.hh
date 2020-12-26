#error "THIS FILE IS NOT BEING USED AND DOES NOT SEEM TO COMPILE!"

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
#ifndef ELLIPT_DGFEMSCHEME_HH
#define ELLIPT_DGFEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/solver/newtoninverseoperator.hh>

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>

#include <dune/fem/schemes/dgelliptic.hh>

template < class Model, int polOrder, SolverType solver >
class DGFemScheme
{
public:
  //! type of the mathematical model
  typedef Model ModelType ;
  typedef typename ModelType::ExactSolutionType ExactSolutionType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R \f$)
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
  // typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Solvers<DiscreteFunctionSpaceType,solver,false> UsedSolverType;
  static_assert( UsedSolverType::solverConfigured, "chosen solver is not configured" );

  typedef typename UsedSolverType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename UsedSolverType::LinearOperatorType LinearOperatorType;

  /*********************************************************/

  //! define Laplace operator
  typedef DifferentiableDGEllipticOperator< LinearOperatorType, ModelType > EllipticOperatorType;

  static const int dimRange = FunctionSpaceType::dimRange;

  DGFemScheme( GridPartType &gridPart,
               const ModelType& implicitModel,
               double penalty,
               const std::string &prefix)
    : implicitModel_( implicitModel ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( prefix.c_str(), discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      implicitOperator_( implicitModel_, discreteSpace_, penalty ),
      linearOperator_( "assempled elliptic operator", discreteSpace_, discreteSpace_ ),
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ),
      penalty_(penalty),
      exactSolution_( implicitModel_.exactSolution(gridPart_) )
  {
    // set all DoF to zero
    solution_.clear();
  }

  DGFemScheme( GridPartType &gridPart,
               const ModelType& implicitModel,
               const std::string &prefix)
    : DGFemScheme( gridPart, implicitModel,
        Dune::Fem::Parameter::getValue<double>("dg.penalty"),
        prefix )
  {
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }
  const ExactSolutionType& exactSolution() const { return exactSolution_; }

  //! solve the system - bool parameter
  //! false: only assemble if grid has changed
  //! true:  assemble in any case
  void solve ( bool assemble )
  {
    typedef typename UsedSolverType::LinearInverseOperatorType LinearInverseOperatorType;
#if 0
    typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, LinearInverseOperatorType > InverseOperatorType;
    InverseOperatorType invOp( implicitOperator_ );
    invOp( rhs_, solution_ );
#else
    implicitOperator_.jacobian( solution_ , linearOperator_ );
    LinearInverseOperatorType invOp( linearOperator_, solverEps_, solverEps_ );
    invOp( rhs_, solution_ );
#endif
  }

protected:
  const ModelType& implicitModel_;   // the mathematical model

  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
  DiscreteFunctionType rhs_;        // the right hand side

  EllipticOperatorType implicitOperator_; // the implicit operator

  LinearOperatorType linearOperator_;  // the linear operator (i.e. jacobian of the implicit)

  const double solverEps_ ; // eps for linear solver
  const double penalty_;
  const ExactSolutionType exactSolution_;
};

#endif // end #if ELLIPT_FEMSCHEME_HH
