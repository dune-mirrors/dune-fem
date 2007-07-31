#include <cmath>
#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;



CG::CG(Communicator &comm) : 
  IterativeLinearSolver(comm), DynamicalObject("CG", comm.id()), 
  r(NULL), d(NULL), h(NULL) 
{}


CG::~CG()
{
  delete[] r;
}


void CG::resize(int new_size, int component)
{
  delete[] r;
  r = new double[3*new_size];
  assert(r);
  d = r + new_size;
  h = r + 2*new_size;
}
 



bool CG::solve(Function &op, double *x, const double *b)
{
  dim = op.dim_of_value();
  new_size(dim);
  int iterations = 0;
  double _tolerance = tolerance;
  if (relative_tolerance){
    double global_dot;
    double local_dot = cblas_ddot(dim, b, 1, b, 1);
    comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);      
    _tolerance *= sqrt(global_dot);
  }

  // preconditioned CG
  if (preconditioner){
    // init
    op(x, r);
    cblas_daxpy(dim, -1.0, b, 1, r, 1);
    (*preconditioner)(r, h);
    cblas_daxpy(dim, -1.0, h, 1, d, 1);
    double local_dot = cblas_ddot(dim, r, 1, h, 1);
    double global_dot;
    comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
    double nu = global_dot;
   
    // iterate
    while (true){
      op(d, h);
      local_dot = cblas_ddot(dim, d, 1, h, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
      const double alpha = nu / global_dot;
      cblas_daxpy(dim, alpha, d, 1, x, 1);
      cblas_daxpy(dim, alpha, h, 1, r, 1);
      (*preconditioner)(r, h);

      iterations++;
      local_dot = cblas_ddot(dim, r, 1, r, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
      if (sqrt(global_dot) < _tolerance 
	  || iterations >= max_num_of_iterations) break;
 
      double beta = 1.0/nu;
      local_dot = cblas_ddot(dim, r, 1, h, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
      nu = global_dot; 
      beta *= nu;
      daxpby(dim, -1.0, h, 1, beta, d, 1);
    }
  }

  // unpreconditioned CG
  else{
    // init
    op(x, r);
    cblas_daxpy(dim, -1.0, b, 1, r, 1);
    cblas_daxpy(dim, -1.0, r, 1, d, 1);
    double local_dot = cblas_ddot(dim, r, 1, r, 1);
    double global_dot;
    comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
    double nu = global_dot;

    // iterate
    while (true){
      op(d, h);
      double local_dot = cblas_ddot(dim, d, 1, h, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
      const double alpha = nu / global_dot;
      cblas_daxpy(dim, alpha, d, 1, x, 1);
      cblas_daxpy(dim, alpha, h, 1, r, 1);

      double beta = 1.0/nu;
      local_dot = cblas_ddot(dim, r, 1, r, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
      nu = global_dot;
      beta *= nu;
      daxpby(dim, -1.0, r, 1, beta, d, 1);
 
      iterations++;
      // std::cout << iterations << " " << sqrt(nu) << std::endl;
      if (sqrt(nu) < _tolerance || iterations >= max_num_of_iterations) break; 
    }
  }
    
 
  // output
  if (IterativeSolver::os){
    *IterativeSolver::os << "CG "<< comm.id() << ": number of iterations: "
			 << iterations
			 << std::endl;
  }
  
  return (iterations < max_num_of_iterations)? true: false;
}




