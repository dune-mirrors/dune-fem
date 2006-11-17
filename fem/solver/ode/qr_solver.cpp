#include <cassert>
#include "linear_solver.hpp"
#include "blas.hpp"


using namespace pardg;



QRSolver::QRSolver() : DynamicalObject("QRSolver", -1), d(NULL) {}


QRSolver::~QRSolver() 
{
  delete[] d;
}


void QRSolver::resize(int new_size, int component)
{
  delete[] d;
  d = new double[new_size]; // new_size >= dim
  assert(d);
}



// QR decomposition of matrix a
// matrix a is provided in row-major format
// first it is converted to column-major format (by transposing)
// to avoid cache mismatches
// possible improvement: tranpose Q at the end of the algorithm
//                       second dot-product in solve routine might
//                       preform better
void QRSolver::prepare(int dim, double *a)
{
  this->a = a;
  this->dim = dim;
  new_size(dim);

  // build transpose of a
  for(int i=0; i<dim; i++){
    for(int j=0; j<i; j++){
      const double tmp = a[i*dim + j];
      a[i*dim + j] = a[j*dim + i];
      a[j*dim + i] = tmp;
    }
  }

  for(int r=0; r<dim-1; r++){
    double *a_rr = a + r*(dim+1);
    const double sigma = cblas_dnrm2(dim-r, a_rr, 1);
    const double beta = 1.0 / sqrt( sigma*( fabs(*a_rr) + sigma ) );

    d[r] = (*a_rr < 0.0)? sigma: -sigma;
    *a_rr -= d[r];
    cblas_dscal(dim-r, beta, a_rr, 1);

    for(int j=r+1; j<dim; j++){
      double *a_rj = a + j*dim + r;
      const double dot = cblas_ddot(dim-r, a_rr, 1, a_rj, 1);
      cblas_daxpy(dim-r, -dot, a_rr, 1, a_rj, 1);
    }
  }

  // last row
  d[dim-1] = a[dim*dim-1];
}



// solve Ax = b <=> QRx = b
// b is initially stored in x
// a: QR decomposition, row-major, d consists of the diagonal elements of R
bool QRSolver::solve(double *x)
{
  // 1. solve Qx = b, with b=x
  for(int r=0; r<dim-1; r++){
    const double dot = cblas_ddot(dim-r, a+r*(dim+1), 1, x+r, 1);
    cblas_daxpy(dim-r, -dot, a+r*(dim+1), 1, x+r, 1);
  }

  // 2. solve Rx = b with b=x
  for(int r=dim-1; r>=0; r--){
    const double dot = cblas_ddot(dim-(r+1), a+r*(dim+1)+dim, dim ,x+(r+1), 1);
    x[r] = (x[r] - dot) / d[r];
  }

  return true;
}
