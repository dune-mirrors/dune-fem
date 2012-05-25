#include <cmath>
#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;

  

BICGSTAB::BICGSTAB(Communicator &comm) :  
  IterativeLinearSolver(comm), DynamicalObject( "BiCGstab", comm.id() ),
  r(NULL), r_star(NULL), p(NULL), s(NULL), tmp(NULL), z(NULL)
{}


BICGSTAB::~BICGSTAB()
{
  delete[] r;
}


void BICGSTAB::set_preconditioner(Function &preconditioner)
{
  this->preconditioner = &preconditioner;
}   


void BICGSTAB::resize(int new_size, int component)
{
  delete[] r;
  r = new double[ ((preconditioner)? 6 : 5) * new_size];
  //r = new double[5 * new_size];
  assert(r);
  r_star = r + new_size;
  p = r + 2*new_size;
  s = r + 3*new_size;
  tmp = r + 4*new_size;
  z = (preconditioner)? r + 5*new_size : 0;
}
 
void BICGSTAB::unset_preconditioner()
{
  this->preconditioner = 0;
}
 
bool BICGSTAB::solve(Function &op, double *x, const double *b)
{
  dim = op.dim_of_value();
  new_size(dim);

  // relative or absolute tolerance
  double local_dot[5], global_dot[5];
  double _tolerance = tolerance;
  if (relative_tolerance){
    local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    _tolerance *= sqrt(global_dot[0]);
  }

  // init
  if (preconditioner){       // right preconditioning
    (*preconditioner)(x, z); // z = M^{-1} x
    op(z, r);                // r = A z
  }
  else op(x, r);             // r = A x

  for(int k=0; k<dim; k++){    
    r[k] = b[k] - r[k];
    p[k] = r[k];
    r_star[k] = r[k];
  }
  local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
  comm.allreduce(1, local_dot, global_dot, MPI_SUM);
  double nu = global_dot[0]; 

 
  // iterate
  int iterations = 0;
  while (true){
    // 2x linear operator 1x dot
    if (preconditioner){    // right preconditioning
      (*preconditioner)(p, z); // z = M^{-1} p
      op(z, tmp);              // tmp = A z
    } 
    else op(p, tmp);           // tmp = A p

    local_dot[0] = cblas_ddot(dim, tmp, 1, r_star, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    const double alpha = nu / global_dot[0];
    for(int k=0; k<dim; k++) s[k] = r[k] - alpha*tmp[k];

    if (preconditioner){       // right preconditioning
      (*preconditioner)(s, z); // z = M^{-1} s
      op(z, r);                // r = A z
    }
    else op(s, r);             // r = A s

    // 5x dot
    local_dot[0]=local_dot[1]=local_dot[2]=local_dot[3]=local_dot[4] = 0.0;
    for(int k=0; k<dim; k++){
      local_dot[0] += r[k]*s[k];
      local_dot[1] += r[k]*r[k];
      local_dot[2] += s[k]*s[k];
      local_dot[3] += s[k]*r_star[k];
      local_dot[4] += r[k]*r_star[k];
    }
    comm.allreduce(5, local_dot, global_dot, MPI_SUM);

    // scalars
    const double omega = global_dot[0] / global_dot[1];
    const double res = sqrt(global_dot[2] 
			    -omega*(2.0*global_dot[0] - omega*global_dot[1]) );

    const double beta = (global_dot[3] - omega*global_dot[4])
      *alpha / (omega*nu);
    
    nu = (global_dot[3] - omega*global_dot[4]);

    // update
    for(int k=0; k<dim; k++){
      x[k] += alpha*p[k] + omega*s[k];
      r[k] = s[k] - omega*r[k];
      p[k] = r[k] + beta*( p[k] - omega*tmp[k] );
    }

    iterations++;    
    if (res < _tolerance || iterations >= max_num_of_iterations) break; 
    if (IterativeSolver::os)
    {
      *IterativeSolver::os << "BiCGstab " << comm.id() << " it: " << iterations << " : " << res << std::endl;
    }
  }
    
 
  // output
  if (IterativeSolver::os){
    *IterativeSolver::os << "BiCGstab " << comm.id() 
			 << ":  number of iterations: " 
			 << iterations
			 << std::endl;
  }

  // setup approx solution for right preconditioner use
  if (preconditioner){ // right preconditioner
    (*preconditioner)(x, z);
    cblas_dcopy(dim, z,1, x,1);
  }

  // update the global number of iterations from IterativeSolver
  num_of_iterations += iterations;

  return (iterations < max_num_of_iterations)? true: false;
}

// bicg withouot preconditioner 
bool BICGSTAB::solve_old(Function &op, double *x, const double *b)
{
  // todo: extend it to preconditioner use
  assert(!preconditioner);

  dim = op.dim_of_value();
  new_size(dim);

  // relative or absolute tolerance
  double local_dot[5], global_dot[5];
  double _tolerance = tolerance;
  if (relative_tolerance){
    local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    _tolerance *= sqrt(global_dot[0]);
  }

  // init
  op(x, r);
  for(int k=0; k<dim; k++){    
    r[k] = b[k] - r[k];
    p[k] = r[k];
    r_star[k] = r[k];
  }
  local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
  comm.allreduce(1, local_dot, global_dot, MPI_SUM);
  double nu = global_dot[0]; 
 
  // iterate
  int iterations = 0;
  while (true){
    // 2x linear operator 1x dot
    op(p, tmp);
    local_dot[0] = cblas_ddot(dim, tmp, 1, r_star, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    const double alpha = nu / global_dot[0];
    for(int k=0; k<dim; k++) s[k] = r[k] - alpha*tmp[k];
    op(s, r);

    // 5x dot
    local_dot[0]=local_dot[1]=local_dot[2]=local_dot[3]=local_dot[4] = 0.0;
    for(int k=0; k<dim; k++){
      local_dot[0] += r[k]*s[k];
      local_dot[1] += r[k]*r[k];
      local_dot[2] += s[k]*s[k];
      local_dot[3] += s[k]*r_star[k];
      local_dot[4] += r[k]*r_star[k];
    }
    comm.allreduce(5, local_dot, global_dot, MPI_SUM);

    // scalars
    const double omega = global_dot[0] / global_dot[1];
    const double res = sqrt(global_dot[2] 
			    -omega*(2.0*global_dot[0] - omega*global_dot[1]) );

    const double beta = (global_dot[3] - omega*global_dot[4])
      *alpha / (omega*nu);
    
    nu = (global_dot[3] - omega*global_dot[4]);

    // update
    for(int k=0; k<dim; k++){
      x[k] += alpha*p[k] + omega*s[k];
      r[k] = s[k] - omega*r[k];
      p[k] = r[k] + beta*( p[k] - omega*tmp[k] );
    }

    iterations++;    
    if (res < _tolerance || iterations >= max_num_of_iterations) break; 
  }
    
 
  // output
  if (IterativeSolver::os){
    *IterativeSolver::os << "BiCGstab " << comm.id() 
			 << ":  number of iterations: " 
			 << iterations
			 << std::endl;
  }

  // update the global number of iterations from IterativeSolver
  num_of_iterations += iterations;

  return (iterations < max_num_of_iterations)? true: false;
}
