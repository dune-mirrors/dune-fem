#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP

#include <iostream>
#include "function.hpp"
#include "iterative_solver.hpp"

namespace DuneODE {

class NonlinearSolver : public IterativeSolver
{
public:
  NonlinearSolver();
  virtual ~NonlinearSolver();

  // solve f(u) = 0,   f:R^n -> R^n some nonlinear function, n=dim
  // result in u
  // return convergence
  virtual bool solve(Function &f, double *u) = 0;

protected:
  int dim;
};



// NonlinearSolver inline implementation
inline
NonlinearSolver::NonlinearSolver() : dim(0) {}

inline
NonlinearSolver::~NonlinearSolver() {}


}

#endif

