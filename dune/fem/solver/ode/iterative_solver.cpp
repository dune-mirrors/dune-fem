#include "iterative_solver.hpp"


IterativeSolver::IterativeSolver() : num_of_iterations(0), os(NULL)
{
  // set this to some useful? values
  tolerance = 1.0e-6;
  relative_tolerance = true;
  max_num_of_iterations = 500;
}



void IterativeSolver::set_tolerance(double tol, bool relative)
{
  tolerance = tol;
  relative_tolerance = relative;
}



void IterativeSolver::set_output(std::ostream &os)
{
  this->os = &os;
}



void IterativeSolver::set_max_number_of_iterations(int iter)
{
  max_num_of_iterations = iter;
}


int IterativeSolver::number_of_iterations() const
{
  return num_of_iterations;
}


void IterativeSolver::reset_number_of_iterations()
{
  num_of_iterations = 0;
}

