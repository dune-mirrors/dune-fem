#include "iterative_solver.hpp"

namespace DuneODE {

IterativeSolver::IterativeSolver() : os(NULL)
{
  // set this to some useful? values
  tolerance = 1.0e-6;
  max_num_of_iterations = 500;
}



void IterativeSolver::set_tolerance(double tol)
{
  tolerance = tol;
}



void IterativeSolver::set_output(std::ostream &os)
{
  this->os = &os;
}



void IterativeSolver::set_max_number_of_iterations(int iter)
{
  max_num_of_iterations = iter;
}

}
