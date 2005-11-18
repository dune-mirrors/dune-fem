#ifndef BLAS_HPP
#define BLAS_HPP


//#define USE_EXTERNAL_BLAS


#ifdef USE_EXTERNAL_BLAS
#include <cblas.h>
#elif USE_INTEL_BLAS
#include <mkl_cblas.h>
#else
#include "subblas.hpp"
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
  while(dim--){
    *x = alpha;
    x += incx;
  } 
}


// y = alpha x + beta y
inline
void daxpby(int dim, 
	    double alpha, const double *x, int incx,
	    double beta, double *y, int incy)
{
  while(dim--){
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
  while(dim--){
    *w = alpha * (*x) + beta * (*y);
    x += incx;
    y += incy;
    w += incw;
  } 
}


#endif
