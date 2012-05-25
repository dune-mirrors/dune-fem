#include <cassert>
#include "ode_solver.hpp"


using namespace pardg;


// class ODESolver 
ODESolver::ODESolver(Communicator &comm, int num_of_tmpobj) : 
  DynamicalObject( "ODESolver", comm.id() ), comm(comm),
  num_of_tmpobj(num_of_tmpobj),
  dim(0), U(NULL), limiter(NULL), os(NULL)
{}


ODESolver::~ODESolver()
{
  delete[] U;
}


void ODESolver::set_output(std::ostream &os)
{
  this->os = &os;
}


void ODESolver::set_limiter(Limiter &limiter)
{
  this->limiter = &limiter;
}


void ODESolver::resize(int new_size, int component)
{
  delete[] U;
  U = new double[new_size*num_of_tmpobj]; // new_size >= dim
  dset(new_size*num_of_tmpobj, 0.0, U, 1);
  assert(U);
}

