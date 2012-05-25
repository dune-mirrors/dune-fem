// FGMRES inline implementation

#include <cmath>
#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;


FGMRES::FGMRES(Communicator &comm, int m) : 
  IterativeLinearSolver(comm), DynamicalObject("FGMRES", comm.id()), 
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



FGMRES::~FGMRES()
{
  delete[] g;
  delete[] v;
  delete[] z;
}


void FGMRES::resize(int new_size, int component)
{
  delete[] v;
  v = new double[(2*m+1)*new_size];
  assert(v);
  z = v + (m+1)*new_size;

  dset((2*m+1)*new_size, 0.0, v, 1);
}



bool FGMRES::solve(Function &op, double *u, const double *b)
{
  // it does not make sense to use FGMRES without a preconditioner
  assert(preconditioner);
  dim = op.dim_of_value();
  new_size(dim);

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
    if (res < tolerance) break;
    g[0] = -res;
    for(int i=1; i<=m; i++) g[i] = 0.0;
    cblas_dscal(dim, 1.0/res, v, 1);

    // iterate
    for(int j=0; j<m; j++){
      double *vj = v + j*dim;
      double *vjp = vj + dim;
      double *zj = z + j*dim;

      // apply preconditioner and the linear operator
      (*preconditioner)(vj, zj);
      op(zj, vjp);      
      
      // Gram-Schmidt procedure
      for(int i=0; i<=j; i++) local_dot[i]=cblas_ddot(dim, vjp, 1, v+i*dim, 1);
      comm.allreduce(j+1, local_dot, global_dot, MPI_SUM);      
      for(int i=0; i<=j; i++) H(i,j) = global_dot[i]; 

      for(int i=0; i<=j; i++) cblas_daxpy(dim, -H(i,j), v+i*dim, 1, vjp, 1);
      local_dot[0] = cblas_ddot(dim, vjp, 1, vjp, 1);
      comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
      H(j+1,j) = sqrt(global_dot[0]); 
      cblas_dscal(dim, 1.0/H(j+1,j), vjp, 1);

      // perform Givens rotation
      for(int i=0; i<j; i++){
        cblas_drot(1, &H(i+1,j), 1, &H(i,j), 1, c[i], s[i]);
      }
      const double h_j_j = H(j,j);
      const double h_jp_j = H(j+1,j);
      const double norm = sqrt(h_j_j*h_j_j + h_jp_j*h_jp_j);
      c[j] = h_j_j / norm;
      s[j] = -h_jp_j / norm;
      cblas_drot(1, &H(j+1,j), 1, &H(j,j), 1, c[j], s[j]);
      cblas_drot(1, &g[j+1], 1, &g[j], 1, c[j], s[j]);

      //*output_stream << fabs(g[j+1]) << std::endl;

      iterations++;
      if (fabs(g[j+1]) < tolerance 
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
      y[i] = (g[i] - dot) / H(i,i);
    }    

    // update the approx. solution
    // u += (z[0], ..., z[last-1]) y
    for(int i=0; i<last; i++) cblas_daxpy(dim, y[i], z+i*dim, 1, u, 1);
    
    if (fabs(g[last]) < tolerance) break;
    if (IterativeSolver::os)
    {
      *IterativeSolver::os<< "FGMRES " << comm.id() 
			    << ": its: " << iterations << "  err: " 
          << fabs(g[last]) << std::endl; 
    }
  }

  // output
  if (IterativeSolver::os){
    *IterativeSolver::os<< "FGMRES " << comm.id() 
			<< ": number of iterations: "
			<< iterations
			<< std::endl;
  }

  return (iterations < max_num_of_iterations)? true: false;
}

 

