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
#ifndef FEM_SOLVER_HH
#define FEM_SOLVER_HH

// C++ includes
#include <iostream>
#include <type_traits>

// include discrete function space
// #include <dune/fem/space/lagrange.hh>
// #include <dune/fem/space/discontinuousgalerkin.hh>

// adaptation ...
// #include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction/managedvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/pardginverseoperators.hh>
#include <dune/fem/solver/oemsolver.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#endif

#if HAVE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#endif

#include <dune/fem/solver/newtoninverseoperator.hh>

// #include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

enum SolverType
{
  matrixFree,  // use the matrix free version of the dune-fem solvers
  fem,         // use the matrix based version of the dune-fem solvers
  femoem,      // use the matrix based version of the dune-fem solvers with blas
  istl,        // use the dune-istl solvers
  umfpack,     // use the direct solver umfpack
  petsc        // use the petsc package
};
enum OperatorType
{
  cg, dg
};

template <class DFSpace, SolverType solver, bool symmetric>
struct Solvers
{
  static const bool solverConfigured = false; // this implementation is used for not installed packages
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // this should work with any discrete function implementation
  typedef Dune::Fem::DynamicVector<double> DofVectorType;
  typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType > > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
};

template <class DFSpace, bool symmetric>
struct Solvers<DFSpace,fem,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename std::conditional<symmetric,
          Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
          Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type
          LinearInverseOperatorType;

};
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,femoem,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // this work with a discrete function implementation based on a double* dof storage
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename std::conditional<symmetric,
          Dune::Fem::OEMCGOp< DiscreteFunctionType, LinearOperatorType >,
          Dune::Fem::OEMBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type
          LinearInverseOperatorType;
};

#if HAVE_DUNE_ISTL
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,istl,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // here we need the special ISTLBlockVectorDiscreteFunction
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename std::conditional<symmetric,
          Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType >,
          Dune::Fem::ISTLGMResOp< DiscreteFunctionType, LinearOperatorType > > :: type
          LinearInverseOperatorType;
};
#endif // HAVE_ISTL
#if HAVE_UMFPACK
template <class DFSpace, bool symmetric>
struct Solvers<DFSpace,umfpack,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType, symmetric > LinearInverseOperatorType;
};
#endif
#if HAVE_PETSC
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,petsc,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
  // to switch between solvers for symmetric and non symmetric operators
  // use the parameter petsc.kspsolver.method
};
#endif

#endif // end #if FEM_SOLVER_HH
