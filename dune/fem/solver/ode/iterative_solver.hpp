#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <iostream>


class IterativeSolver
{
public:
  void set_tolerance(double tol, bool relative = true);
  void set_max_number_of_iterations(int iter);
  void set_output(std::ostream &os);
  int number_of_iterations() const;
  void reset_number_of_iterations();

protected:
  IterativeSolver();

  double tolerance;  
  bool relative_tolerance;
  int max_num_of_iterations, num_of_iterations;  
  std::ostream *os;
};



#endif
