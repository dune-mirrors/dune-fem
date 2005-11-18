#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <iostream>

namespace DuneODE {

class IterativeSolver
{
public:
  void set_tolerance(double tol);
  void set_max_number_of_iterations(int iter);
  void set_output(std::ostream &os);

protected:
  IterativeSolver();

  double tolerance;
  int max_num_of_iterations;
  std::ostream *os;
};

}


#endif
