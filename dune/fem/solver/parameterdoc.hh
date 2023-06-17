#ifndef DUNE_FEM_SOLVER_PARAMETERDOC_HH
#define DUNE_FEM_SOLVER_PARAMETERDOC_HH

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <set>
#include <vector>
#include <iomanip>

#include <dune/fem/solver/parameter.hh>

namespace Dune
{
  namespace Fem{

  namespace detail {

        std::pair< std::string, std::string > solverString(const bool havePetsc);
  } // namespace

  }
}

template <class Scheme>
auto docString()
{
  bool havePetsc = false;
#if HAVE_PETSC
  havePetsc = true;
#endif
  auto str = Dune::Fem::detail::solverString(havePetsc);
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
