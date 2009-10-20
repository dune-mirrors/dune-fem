#ifndef BLAS_HPP
#define BLAS_HPP

#include <cstring>
#ifdef USE_EXTERNAL_BLAS
// use atlas cblas 
extern "C" {
#include <atlas/cblas.h>
}
#else
// use own blas implementation
#include <dune/fem/solver/ode/subblas.hpp>
#endif

// non standard functions (?)

// y = alpha x + beta y
void daxpby(int dim, 
	    double alpha, const double *x, int incx,
	    double beta, double *y, int incy);

// w = alpha x + beta y
void dwaxpby(int dim, 
	     double alpha, const double *x, int incx,
	     double beta, const double *y, int incy, 
	     double *w, int incw);


// extra functions

// x[i] = alpha, i=0,..,dim-1
void dset(int dim, double alpha, double *x, int incx);





// ========= inline implementation

// x[i] = alpha, i=0,..,dim-1
inline
void dset(int dim, double alpha, double *x, int incx)
{
  if (incx == 1 && alpha == 0.0){
    memset(x, 0, dim*sizeof(double));
  }
  else{
    int i = dim;
    while(i--){
      *x = alpha;
      x += incx;
    } 
  }
}


// y = alpha x + beta y
inline
void daxpby(int dim, 
	    double alpha, const double *x, int incx,
	    double beta, double *y, int incy)
{
  int i = dim;
  while(i--){
    *y = alpha * (*x) + beta * (*y);
    x += incx;
    y += incy;
  } 
}


// w = alpha x + beta y
inline
void dwaxpby(int dim, 
	     double alpha, const double *x, int incx,
	     double beta, const double *y, int incy, 
	     double *w, int incw)
{
  int i = dim;

  while(i--){
    *w = alpha * (*x) + beta * (*y);
    x += incx;
    y += incy;
    w += incw;
  } 
}


#endif
