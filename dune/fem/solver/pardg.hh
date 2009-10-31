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

//#ifdef ENABLE_PARDG 
#define USE_PARDG_ODE_SOLVER 
//#endif

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
#ifdef ENABLE_PARDG 
#include <timer.hpp>
#else 
#include "ode/timer.hpp"
#endif
} // end namespace pardg 

// include pardg communicator 
#include "ode/communicator.hpp"
} // end namespace PARDG_NS

// if pardg library was found 
#ifdef ENABLE_PARDG 

#include <blas.hpp>
namespace PARDG_NS {
namespace pardg {
// we also need vector to be in namespace parDG 
#include <vector.hpp>
} // end namespace pardg 

#include <quadrature.hpp>  
#include <ode_solver.hpp>
#include <linear_solver.hpp>
} // end namespace PARDG_NS
// else use build in ode solver (may be outdated)
#else

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
} // end namespace PARDG_NS
#endif

#endif // end USE_DENNIS 

#endif
