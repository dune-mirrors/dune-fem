#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <iostream>
#include "dynamical_object.hpp"
#include "function.hpp"
#include "nonlinear_solver.hpp"


class Newton : public NonlinearSolver, public DynamicalObject
{
public:
  Newton();
  ~Newton();

  // from NonlinearSolver
  virtual bool solve(Function &F, double *u);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  double *f, *Df;
  int *p; // permutation array for LU decomp
};


#endif
