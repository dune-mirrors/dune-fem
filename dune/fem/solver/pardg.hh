#ifndef PARDG_INCLUDE_HH
#define PARDG_INCLUDE_HH

#include <unistd.h>
#include <time.h>
#include <sys/times.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

//#ifdef ENABLE_PARDG 
#define USE_PARDG_ODE_SOLVER 
//#endif

// if the preprocessor variable is defined, the ODE Solver from Dennis
// are used.
#ifdef USE_PARDG_ODE_SOLVER

// timer has no namespace therefore we put here 
namespace pardg {
// if pardg library was found 
#ifdef ENABLE_PARDG 
#include <timer.hpp>
#else 
#include "ode/timer.hpp"
#endif
}

// include pardg communicator 
#include "ode/communicator.hpp"

// if pardg library was found 
#ifdef ENABLE_PARDG 

#include <blas.hpp>
namespace pardg {
// we also need vector to be in namespace parDG 
#include <vector.hpp>
}

#include <quadrature.hpp>  
#include <ode_solver.hpp>
#include <linear_solver.hpp>

// else use build in ode solver (may be outdated)
#else

namespace pardg {
// we also need vector to be in namespace parDG 
#include "ode/vector.hpp"
}

#include "ode/quadrature.hpp"  
#include "ode/blas.hpp"
#include "ode/function.hpp"
#include "ode/ode_solver.hpp"
#include "ode/linear_solver.hpp"
#endif

#endif // end USE_DENNIS 

#endif
