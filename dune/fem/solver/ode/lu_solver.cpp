#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;



LUSolver::LUSolver() : DynamicalObject("LUSolver", -1), p(NULL) {}


LUSolver::~LUSolver() 
{
  delete[] p;
}


void LUSolver::resize(int new_size, int component)
{
  delete[] p;
  p = new int[new_size]; // new_size >= dim
  assert(p);
}



// LU decomposition of matrix a
void LUSolver::prepare(int dim, double *a)
{
  this->a = a;
  this->dim = dim;
  new_size(dim);

  for(int i=0; i<dim-1; i++){
    // pivot search
    p[i] = cblas_idamax(dim-i, a+i*(dim+1), dim) + i;

    // exchange row i with row argmax=p[i]
    cblas_dswap(dim, a+i*dim, 1, a+p[i]*dim, 1);

    // elimination
    const double *a_i = a + i*(dim+1);
    for(int k=i+1; k<dim; k++){
      double *a_k = a + k*dim + i;
      const double lambda = *a_k / *a_i;
      *a_k = lambda;
      cblas_daxpy(dim-(i+1), -lambda, a_i+1, 1, a_k+1, 1);
    }
  }

  // last row
  p[dim-1] = dim-1;
}



// solve Ax = b <=> LUx = Pb
// b is initially stored in x
// a: LU decomposition, row-major format, p is the permutation array
bool LUSolver::solve(double *x)
{
  // 1. x = Px, permutation with right hand side
  for(int i=0; i<dim-1; i++){
    double tmp = x[i];
    x[i] = x[p[i]];
    x[p[i]] = tmp;
  }

  // 2. Lx = x, forward solve
  for(int i=0; i<dim; i++) x[i] -= cblas_ddot(i, a+i*dim, 1, x, 1);

  // 3. Ux = x, backward solve
  for(int i=dim-1; i>=0; i--){
    const double *a_ii = a+i*(dim+1);
    x[i] = ( x[i] - cblas_ddot(dim-(i+1), a_ii+1, 1, x+i+1, 1) ) / (*a_ii);
  }

  return true;
}
