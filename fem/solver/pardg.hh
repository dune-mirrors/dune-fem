#ifndef PARDG_INCLUDE_HH
#define PARDG_INCLUDE_HH

#define USE_PARDG_ODE_SOLVER

// if the preprocessor variable is defined, the ODE Solver from Dennis
// are used.
#ifdef USE_PARDG_ODE_SOLVER

// timer has no namespace therefore we put here 
namespace pardg {
#include "ode/timer.hpp"
}

// include pardg communicator 
#include "ode/communicator.hpp"

// if pardg library was found 
#ifdef ENABLE_PARDG 

#include <ode_solver.hpp>
#include <linear_solver.hpp>

// else use build in ode solver (may be outdated)
#else

#include "ode/blas.hpp"
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"

#endif

#endif // end USE_DENNIS 

#endif
