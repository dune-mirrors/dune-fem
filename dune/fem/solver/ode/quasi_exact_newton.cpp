#include <cmath> 
#include <cassert>
#include <cfloat>
#include "quasi_exact_newton.hpp"
#include "blas.hpp"

using namespace pardg;



QuasiExactNewton::QuasiExactNewton(Communicator &comm) : 
  DynamicalObject( "QuasiExactNewton", comm.id() ),
  f(NULL), p(NULL), u_tmp(NULL), op(*this), comm(comm), linear_solver(NULL) 
{
  // set this to some useful values
  tolerance = 1.0e-6;
  max_num_of_iterations = 20;
}


QuasiExactNewton::~QuasiExactNewton()
{
  delete[] f;
}


void QuasiExactNewton::resize(int new_size, int component)
{
  delete[] f;
  f = new double[3 * new_size]; // new_size >= dim
  assert(f);
  p = f + new_size;
  u_tmp = p + new_size;
}


void QuasiExactNewton::set_linear_solver(IterativeLinearSolver &ls)
{
  linear_solver = &ls;
}


bool QuasiExactNewton::solve(Function &F, double *u)
{
  dim = F.dim_of_argument();
  new_size(dim);
  this->u = u;
  this->F = &F;
  int num_of_iterations = 0;
  double res_old = DBL_MAX;

  while (num_of_iterations < max_num_of_iterations) {
    F(u, f);
    
    double local_dot, global_dot;
    local_dot = cblas_ddot(dim, f, 1, f, 1);
    comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
    const double res = sqrt(global_dot);
    if (res < tolerance) return true;
    if (res >= res_old) return false;
    res_old = res;

    assert(linear_solver);
    dset(dim, 0.0, p, 1);
    bool lin_solver_conv = linear_solver->solve(op, p, f);
    if (!lin_solver_conv) return false;

    cblas_daxpy(dim, -1.0, p, 1, u, 1);

    if (NonlinearSolver::os){    
      *NonlinearSolver::os << "QuasiExactNewton: iteration: " 
			   << num_of_iterations << "    "
			   << "res: " << res << "   "
			   << std::endl;
    }

    num_of_iterations++;
  } 
  
  return false; // num_of_iterations >= max_num_of_iterations
}






// ======== QuasiExactNewton::LinearOperator implementation
QuasiExactNewton::LinearOperator::LinearOperator(QuasiExactNewton &qen) :
  qen(qen) {}


int QuasiExactNewton::LinearOperator::dim_of_argument(int i) const
{
  return qen.dim;
}


int QuasiExactNewton::LinearOperator::dim_of_value(int i) const
{
  return qen.dim;
}


void QuasiExactNewton::LinearOperator::operator()(const double *p, 
						  double *DFu_p, int i)
{
  double local_dot[2], global_dot[2];
  local_dot[0] = cblas_ddot(qen.dim, qen.u, 1, qen.u, 1);
  local_dot[1] = cblas_ddot(qen.dim, p, 1, p, 1);
  qen.comm.allreduce(2, local_dot, global_dot, MPI_SUM);
  double norm_u = sqrt(global_dot[0]);
  double norm_p_sq = global_dot[1];

  double eps = (norm_p_sq > DBL_EPSILON)? 
    sqrt( (1.0+norm_u)*DBL_EPSILON / norm_p_sq ) : sqrt(DBL_EPSILON);
  
  //for(int k=0; k<qen.dim; k++) qen.u_tmp[k] = qen.u[k] + eps*p[k];
  dwaxpby(qen.dim, 1.0, qen.u, 1, eps, p, 1, qen.u_tmp, 1);
  (*qen.F)(qen.u_tmp, DFu_p);
  //for(int k=0; k<qen.dim; k++) DFu_p[k] = (DFu_p[k] - qen.f[k])/eps;
  daxpby(qen.dim, -1.0/eps, qen.f, 1, 1.0/eps, DFu_p, 1);
}
