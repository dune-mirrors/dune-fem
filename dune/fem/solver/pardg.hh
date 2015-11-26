#ifndef PARDG_INCLUDE_HH
#define PARDG_INCLUDE_HH

#include <unistd.h>
#include <time.h>
#include <sys/times.h>

#include <cassert>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <iostream>

#define USE_PARDG_ODE_SOLVER

// use different namespaces in case of MPI or not
#if HAVE_MPI
#include <mpi.h>
#define PARDG_NS parDG_MPI
#else
#define PARDG_NS parDG_NoMPI
#endif

// define combined namespace for later use
#define PARDG PARDG_NS::pardg

// if the preprocessor variable is defined, the ODE Solver from Dennis
// are used.
#ifdef USE_PARDG_ODE_SOLVER

// timer has no namespace therefore we put here
namespace PARDG_NS {
namespace pardg {
// if pardg library was found
#include "ode/timer.hpp"
} // end namespace pardg

// include pardg communicator
#include "ode/communicator.hpp"
} // end namespace PARDG_NS

namespace PARDG_NS {
namespace pardg {
// we also need vector to be in namespace parDG
#include "ode/vector.hpp"
} // end namespace pardg
#include "ode/blas.hpp"
#include "ode/quadrature.hpp"
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"
  void set_tolerance(IterativeSolver &solver,
      double redEps, double absLimit, std::string &paramName)
  {
    static const std::string errorTypeTable[] =
      { "absolute", "relative", "residuumReduction" };
    int errorType = Parameter::getEnum( paramName, errorTypeTable, 0 );
    switch (errorType)
    {
      case 1: solver.set_tolerance(absLimit, ToleranceCriteria::absolute); break;
      case 2: solver.set_tolerance(redEps, ToleranceCriteria::relative); break;
      case 3: solver.set_tolerance(redEps, ToleranceCriteria::residuumReduction); break;
    }
  }
  void set_tolerance(IterativeSolver &solver, double tol, std::string &paramName)
  {
    set_tolerance(solver,tol,tol,paramName);
  }
} // end namespace PARDG_NS

#endif // end USE_PARDG_ODE_SOLVER

#endif // #ifndef PARDG_INCLUDE_HH
