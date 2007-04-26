// GMRES inline implementation

#ifndef ODE_GMRES_CPP 
#define ODE_GMRES_CPP 

#include <cmath>
#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"

using namespace pardg; 

GMRES::GMRES(Communicator &comm, int m) : 
  IterativeLinearSolver(comm), DynamicalObject("GMRES", comm.id()), 
  m(m), H(m+1,m), v(NULL), z(NULL)
{
  g = new double[6*m+1];
  assert(g);
  s = g + (m+1);
  c = s + m;
  y = c + m;
  local_dot = y + m;
  global_dot = local_dot + m;

  dset(6*m+1, 0.0, g, 1);
}



GMRES::~GMRES()
{
  delete[] g;
  delete[] v;
  delete[] z;
}



void GMRES::set_preconditioner(Function &preconditioner)
{
  unset_preconditioner();

  this->preconditioner = &preconditioner;
  if (dim > 0){
    z = new double[size()];
    assert(z);
  }
}

void GMRES::unset_preconditioner()
{
  this->preconditioner = 0; 
  delete [] z;
  z = 0;
}


void GMRES::resize(int new_size, int component)
{
  if (preconditioner){
    delete[] z;
    z = new double[new_size];
    assert(z);
    dset(new_size, 0.0, z, 1);
  }

  delete[] v;
  v = new double[(m+1)*new_size];
  assert(v);

  dset((m+1)*new_size, 0.0, v, 1);
}
 



bool GMRES::solve(Function &op, double *u, const double *b)
{
  dim = op.dim_of_value();
  new_size(dim);

  // relative or absolute tolerance
  double _tolerance = tolerance;
  if (relative_tolerance){
    local_dot[0] = cblas_ddot(dim, b, 1, b, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    _tolerance *= sqrt(global_dot[0]);
  }

  int iterations = 0;
  while (true){
    // start
    op(u, v);
    cblas_daxpy(dim, -1.0, b, 1, v, 1);
    local_dot[0] = cblas_ddot(dim, v, 1, v, 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    double res = sqrt(global_dot[0]); 
    if (res < _tolerance) break;
    g[0] = -res;
    for(int i=1; i<=m; i++) g[i] = 0.0;
    cblas_dscal(dim, 1.0/res, v, 1);
    
    // iterate
    for(int j=0; j<m; j++){
      double *vj = v + j*dim;
      double *vjp = vj + dim;

      // apply the linear operator (perhaps in combination with the 
      // preconditioner)
      if (preconditioner)
      {
	      (*preconditioner)(vj, z);
      	op(z, vjp);
      }
      else op(vj, vjp);

      cblas_dgemv(CblasRowMajor, CblasNoTrans, 
		  j+1, dim, 1.0, v, dim, vjp, 1, 0.0, local_dot, 1);
 
      comm.allreduce(j+1, local_dot, global_dot, MPI_SUM);      
      for(int i=0; i<=j; i++) H(i,j) = global_dot[i]; 



      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  1, dim, j+1,  -1.0, global_dot, m,  v, dim,  1.0, vjp, dim);

      // should be replaced by the following, but doesnt work
      // ???
//       cblas_dgemv(CblasRowMajor, CblasTrans,
// 		  dim, j+1, -1.0, v, dim,  global_dot, 1,  1.0, vjp, 1);

      local_dot[0] = cblas_ddot(dim, vjp, 1, vjp, 1);
      comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
      H(j+1,j) = sqrt(global_dot[0]); 
      cblas_dscal(dim, 1.0/H(j+1,j), vjp, 1);

      // perform Givens rotation
      for(int i=0; i<j; i++)
      {
      	cblas_drot(1, &H(i+1,j), 1, &H(i,j), 1, c[i], s[i]);
      }
      const double h_j_j = H(j,j);
      const double h_jp_j = H(j+1,j);
      const double norm = sqrt(h_j_j*h_j_j + h_jp_j*h_jp_j);
      c[j] = h_j_j / norm;
      s[j] = -h_jp_j / norm;
      cblas_drot(1, &H(j+1,j), 1, &H(j,j), 1, c[j], s[j]);
      cblas_drot(1, &g[j+1], 1, &g[j], 1, c[j], s[j]);

      //*os << fabs(g[j+1]) << std::endl;

      iterations++;
      if (fabs(g[j+1]) < _tolerance 
	        || iterations >= max_num_of_iterations) break;
    }
    
    //
    // form the approximate solution
    //

    int last = iterations%m;
    if (last == 0) last = m;
    // compute y via backsubstitution
    for(int i=last-1; i>=0; i--){
      const double dot = cblas_ddot(last-(i+1), &H(i,i)+1, 1, &y[i+1], 1);
      y[i] = (g[i] - dot)/ H(i,i);
    }    

    // update the approx. solution
    if (preconditioner)
    {
      // u += M^{-1} (v[0], ..., v[last-1]) y	
      double *u_tmp = v + m*dim; // we don't need this vector anymore
      dset(dim, 0.0, u_tmp, 1);
      for(int i=0; i<last; i++)
      {
	      double *vi = v + i*dim;
	      cblas_daxpy(dim, y[i], vi, 1, u_tmp, 1);
      }
      (*preconditioner)(u_tmp, z);
      cblas_daxpy(dim, 1.0, z, 1, u, 1);
    }
    else{
      // u += (v[0], ..., v[last-1]) y
      for(int i=0; i<last; i++)
      {
      	double *vi = v + i*dim;
      	cblas_daxpy(dim, y[i], vi, 1, u, 1);
      }
    }

    if (fabs(g[last]) < _tolerance) break;
    if (IterativeSolver::os)
    {
      *IterativeSolver::os << "GMRES "<< comm.id() << " it: " << iterations << " : " <<  fabs(g[last]) << std::endl; 
    }

  }

  // output
  if (IterativeSolver::os){
    *IterativeSolver::os << "GMRES " << comm.id() 
			 << ": number of iterations: "      
			 << iterations
			 << std::endl;
  }

  // update the global number of iterations from IterativeSolver
  num_of_iterations += iterations;

  return (iterations < max_num_of_iterations)? true: false;
}

#endif
