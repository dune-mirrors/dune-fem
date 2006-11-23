#include <cmath>
#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;

  

BICGSTAB::BICGSTAB(Communicator &comm) :  
  IterativeLinearSolver(comm), DynamicalObject( "BiCGstab", comm.id() ),
  r(NULL), r_star(NULL), p(NULL), s(NULL), tmp(NULL) , z(NULL)
{}


BICGSTAB::~BICGSTAB()
{
  delete[] r;
  delete[] z;
}

void BICGSTAB::set_preconditioner(Function &preconditioner)
{
  this->preconditioner = &preconditioner;
  if (dim > 0)
  {
    delete [] z;
    z = new double[size()];
    assert(z);
    dset(size(), 0.0, z, 1);
  }
}

void BICGSTAB::unset_preconditioner()
{
  this->preconditioner = 0; 
  if( z ) delete [] z;  
  z = 0;
}

void BICGSTAB::resize(int new_size, int component)
{
  delete[] r;
  r = new double[5*new_size];
  assert(r);
  r_star = r + new_size;
  p = r + 2*new_size;
  s = r + 3*new_size;
  tmp = r + 4*new_size;

  if(this->preconditioner) 
  {
    delete [] z;
    z = new double [new_size];
    dset(new_size, 0.0, z, 1);
  }
}
 

bool BICGSTAB::solve_old(Function &op, double *x, const double *b)
{
  // todo: extend it to preconditioner use
  assert(!preconditioner);

  dim = op.dim_of_value();
  new_size(dim);

  // init
  op(x, r);
  for(int k=0; k<dim; k++){    
    r[k] = b[k] - r[k];
    p[k] = r[k];
    r_star[k] = r[k];
  }
  double local_dot[2], global_dot[2];
  local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
  comm.allreduce(1, local_dot, global_dot, MPI_SUM);
  double nu = global_dot[0]; 
 
  // iterate
  int iterations = 0;
  while (true){
    op(p, tmp);
    local_dot[0] = cblas_ddot(dim, tmp, 1, r_star, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    const double alpha = nu / global_dot[0];
    for(int k=0; k<dim; k++) s[k] = r[k] - alpha*tmp[k];
    op(s, r);

    // dot[0]=(s,r)=(s,As), dot[1]=(r,r)=(As,As)
    local_dot[0] = local_dot[1] = 0.0;
    for(int k=0; k<dim; k++){
      local_dot[0] += r[k]*s[k];
      local_dot[1] += r[k]*r[k];
    }
    comm.allreduce(2, local_dot, global_dot, MPI_SUM);
    const double omega = global_dot[0] / global_dot[1];

    for(int k=0; k<dim; k++){
      x[k] += alpha*p[k] + omega*s[k];
      r[k] = s[k] - omega*r[k];
    }

    local_dot[0] = cblas_ddot(dim, r, 1, r, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    double res = sqrt(global_dot[0]);

    if (IterativeSolver::os && comm.id() == 0){
      *IterativeSolver::os << "BiCGstab " << comm.id() << ": iteration: " 
		     << iterations << "     " 
		     << "res: " << res 
		     << std::endl;
    }
    iterations++;    
    if (res < tolerance || iterations >= max_num_of_iterations) break; 

    double beta = alpha / (omega * nu);
    local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    nu = global_dot[0];
    beta *= nu;

    for(int k=0; k<dim; k++) p[k] = r[k] + beta*( p[k] - omega*tmp[k] );
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




bool BICGSTAB::solve(Function &op, double *x, const double *b)
{
  // todo: extend it to preconditioner use
  //assert(!preconditioner);

  dim = op.dim_of_value();
  new_size(dim);

  // relative or absolute tolerance
  double local_dot[5], global_dot[5];
  double _tolerance = tolerance;

  if (relative_tolerance){
    //local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    if (this->preconditioner)
    {
      (*(this->preconditioner))(b, z);
      local_dot[0] = cblas_ddot(dim, z, 1, z, 1);
    }
    else 
    {
      local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    }
    //local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    _tolerance *= sqrt(global_dot[0]);
  }

  // init
  if (this->preconditioner)
  {
    (*(this->preconditioner))(x, z);
    op(z, r);
  }
  else 
  {
    op(x, r);
  }

  if(this->preconditioner)
  {
    (*(this->preconditioner))(b, z); 
  }
  
  const double * rhs = (this->preconditioner) ? z : b;
  for(int k=0; k<dim; k++)
  {    
    r[k] = rhs[k] - r[k];
    p[k] = r[k];
    r_star[k] = r[k];
  }
  
  local_dot[0] = cblas_ddot(dim, r, 1, r_star, 1);
  comm.allreduce(1, local_dot, global_dot, MPI_SUM);
  double nu = global_dot[0]; 
 
  // iterate
  int iterations = 0;
  while (true)
  {
    // 2x linear operator 1x dot
    if (this->preconditioner)
    {
      (*(this->preconditioner))(p, z);
      op(z, tmp);
    }
    else 
    {
      op(p, tmp);
    }

    local_dot[0] = cblas_ddot(dim, tmp, 1, r_star, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);
    const double alpha = nu / global_dot[0];
    for(int k=0; k<dim; k++) s[k] = r[k] - alpha*tmp[k];
    
    if (this->preconditioner)
    {
      (*(this->preconditioner))(s, z);
      op(z, r);
    }
    else 
    {
      op(s, r);
    }

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




