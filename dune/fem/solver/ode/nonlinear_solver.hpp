#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP

#include <iostream>
#include "function.hpp"
#include "iterative_solver.hpp"


namespace pardg
{


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


} // namespace pardg




// NonlinearSolver inline implementation
inline
pardg::NonlinearSolver::NonlinearSolver() : dim(0) {}

inline
pardg::NonlinearSolver::~NonlinearSolver() {}




#endif

