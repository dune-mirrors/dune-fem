#ifndef DUNE_FEM_SOLVER_PARAMETERDOC_HH
#define DUNE_FEM_SOLVER_PARAMETERDOC_HH

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <set>
#include <vector>
#include <iomanip>

#include <dune/fem/operator/matrix/colcompspmatrix.hh>
#include <dune/fem/operator/matrix/istlpreconditioner.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/solver/petscinverseoperators.hh>
#include <dune/fem/solver/parameter.hh>

namespace {

      template <typename function_t>
      std::set< std::string > addMethods( const std::vector<int>& mth, function_t* fct )
      {
        std::set< std::string > s;
        for( const auto& m : mth )
          s.insert( fct( m ) );
        return s;
      }

      template <class DiscreteFunction>
      std::pair< std::string, std::string > solverString()
      {
        std::set< std::string > numpy, istl, petsc, all;
        std::set< std::string > preNum, preIstl, prePetsc, preAll;

        using Dune::Fem::SolverParameter;

        // numpy
        {
          typedef Dune::Fem::KrylovInverseOperator< DiscreteFunction >  InverseOperatorType;
          numpy = addMethods( InverseOperatorType::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          preNum = addMethods( InverseOperatorType::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          for( const auto& m : numpy )
            all.insert( m );
          for( const auto& m : preNum )
            preAll.insert( m );
        }

#if HAVE_DUNE_ISTL
        // istl
        {
          istl = addMethods( Dune::Fem::ISTLInverseOperatorMethods::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          preIstl = addMethods( Dune::Fem::ISTLPreconditionMethods::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          for( const auto& m : istl )
            all.insert( m );
          for( const auto& m : preIstl )
            preAll.insert( m );
        }
#endif

#if HAVE_PETSC
        // petsc
        {
          typedef Dune::Fem::PetscInverseOperator< DiscreteFunction > InverseOperatorType;
          petsc = addMethods( InverseOperatorType::supportedSolverMethods(), &SolverParameter::solverMethodTable );
          for( const auto& m : petsc )
            all.insert( m );
          prePetsc = addMethods( InverseOperatorType::supportedPreconditionMethods(), &SolverParameter::preconditionMethodTable );
          auto extra = InverseOperatorType::extraPreconditionMethods();
          for( const auto& p : extra )
            prePetsc.insert( p );
          for( const auto& m : prePetsc )
            preAll.insert( m );
        }
#endif

        const auto contains = [&](const std::string& m, const std::set< std::string >& s ) -> std::string
        {
          std::string yes(" x ");
          if( m == "ilu" || m == "lu" )
            yes = std::string(" s ");
          auto it = s.find( m );
          if( it == s.end() )
            return std::string("---");
          else
            return yes;
        };

        std::stringstream out;
        out << "------------------------------------------" << std::endl;
        out << "|  Solver  |         Storage             |" << std::endl;
        out << "|   name   |  numpy  |  istl   |  petsc  |" << std::endl;
        out << "|----------|---------|---------|---------|" << std::endl;
        for( const auto & m : all )
        {
          out << "| " << std::setw(8) << std::left << m << " |   " << contains(m,numpy)
                                                        << "   |   " << contains(m,istl)
                                                        << "   |   " << contains(m,petsc) << "   |" << std::endl;
        }
        out << "------------------------------------------" << std::endl;

        std::stringstream pre;
        pre << "-----------------------------------------------" << std::endl;
        pre << "|  Precondition | (x = parallel | s = serial) |" << std::endl;
        pre << "|  method       |  numpy  |   istl  |  petsc  |" << std::endl;
        pre << "|---------------|---------|---------|---------|" << std::endl;
        for( const auto & m : preAll )
        {
          pre << "| " << std::setw(12) << std::left << m << "  |   " << contains(m,preNum)
                                                        << "   |   " << contains(m,preIstl)
                                                        << "   |   " << contains(m,prePetsc) << "   |" << std::endl;
        }
        pre << "-----------------------------------------------" << std::endl;

        return std::make_pair( out.str(), pre.str() );
      }
} // namespace


template <class Scheme>
auto docString()
{
  auto str = solverString< typename Scheme::DiscreteFunctionType >();
  std::string doc = R"doc(
A scheme finds a solution `u=ufl.TrialFunction` for a given variational equation.
The main method is `solve` which takes a discrete functions as `target` argument to
store the solution. The method always uses a Newton method to solve the problem.
The linear solver used in each iteration of the Newton method can be chosen
using the `solver` parameter in the constructor of the scheme. Available solvers are:
)doc"
+ str.first +
R"doc(
In addition the direct solvers from the `suitesparse` package can be used with the
`numpy` storage. In this case provide a tuple as `solver` argument with "suitesparse" as
first argument and the solver to use as second, e.g.,
'solver=("suitesparse","umfpack")'.

The detailed behavior of the schemes can be customized by providing a
`parameters` dictionary to the scheme constructor, e.g.,
   {"newton.tolerance": 1e-3, # tolerance for newton solver
    "newton.verbose": False,  # toggle iteration output
    "newton.linear.tolerance": 1e-5, # tolerance for linear solver
    "newton.linear.errormeasure": "absolute", # or "relative" or "residualreduction"
    "newton.linear.preconditioning.method": "jacobi", # (see table below)
    "newton.linear.preconditioning.hypre.method": "boomeramg", #  "pilu-t" "parasails"
    "newton.linear.preconditioning.iteration": 3, # iterations for preconditioner
    "newton.linear.preconditioning.relaxation": 1.0, # omega for SOR and ILU
    "newton.linear.maxiterations":1000, # max number of linear iterations
    "newton.linear.verbose": False,     # toggle linear iteration output
    "newton.linear.preconditioning.level": 0} # fill-in level for ILU preconditioning
)doc"
+ str.second +
R"doc(
The functionality of some of the preconditioners listed for petsc will
depend on the petsc installation.
)doc";
return doc;
 }
#endif
