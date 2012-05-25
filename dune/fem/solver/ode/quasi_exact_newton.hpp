#ifndef QUASI_EXACT_NEWTON_HPP
#define QUASI_EXACT_NEWTON_HPP

#include <iostream>
#include "communicator.hpp"
#include "dynamical_object.hpp"
#include "function.hpp"
#include "nonlinear_solver.hpp"
#include "linear_solver.hpp"


namespace pardg
{


class QuasiExactNewton : public NonlinearSolver, public DynamicalObject
{
public:
  class LinearOperator : public Function
  {
  public:
    LinearOperator(QuasiExactNewton &qen);

    // from Function
    virtual void operator()(const double *p, double *DFu_p, int i);
    virtual int dim_of_argument(int i) const;
    virtual int dim_of_value(int i) const;

  private:
    QuasiExactNewton &qen;
  };


  QuasiExactNewton(Communicator &comm);
  ~QuasiExactNewton();
  void set_linear_solver(IterativeLinearSolver &ls);

  // from NonlinearSolver
  virtual bool solve(Function &F, double *u);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  friend class LinearOperator;

  double *f, *p, *u_tmp; // u_tmp provides a temp-vector for LinearOperator
  LinearOperator op;

  double *u;
  Function *F;

  Communicator &comm;
  IterativeLinearSolver *linear_solver;
};


} // namespace pardg


#endif
