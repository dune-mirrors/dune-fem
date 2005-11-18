#ifndef SUBBLAS_HPP
#define SUBBLAS_HPP

#include <cmath>
#include <cassert>


// a subset of BLAS subroutines


// Enumerated and derived types
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};


// BLAS level 1 functions
void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drot(int dim, double *x, int incx, double *y, int incy, 
		double c, double s);
void cblas_dswap(int dim, double *x, int incx, double *y, int incy);
void cblas_dcopy(int dim, const double *x, int incx, double *y, int incy);
void cblas_dscal(int dim, double alpha, double *x, int incx);
void cblas_daxpy(int dim, double alpha, const double *x, int ix, 
		 double *y, int iy);
double cblas_ddot(int dim, const double *x, int incx, 
		  const double *y, int incy);
double cblas_dnrm2(int dim, const double *x, int incx);
int cblas_idamax(int dim, const double *x, int incx);


// BLAS level 2 functions
void cblas_dgemv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans,
		 int m, int n, double alpha, const double *A, int lda,
		 const double *x, int incx, double beta, double *y, int incy);






// ========= inline implementation

// some BLAS Level 1 functions -- inline implementation
inline
void cblas_drotg(double *sa, double *sb, double *c, double *s)
{
  double roe = *sb;
  if( fabs(*sa) > fabs(*sb) ) roe = *sa;
  double scale = fabs(*sa) + fabs(*sb);
  double r, z;

  if( scale == 0.0 ){
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
    z = 0.0;
  }
  else{
    r = scale * sqrt( (*sa/scale)*(*sa/scale) + (*sb/scale)*(*sb/scale) );
    if (roe < 0.0)  r = -r;
    *c = *sa/r;
    *s = *sb/r;
    z = 1.0;
    if( fabs(*sa) > fabs(*sb) ) z = *s;
    if( fabs(*sb) >= fabs(*sa) && *c != 0.0 ) z = 1.0/(*c);
  }

  *sa = r;
  *sb = z;
}


inline
void cblas_drot(int dim, double *x, int incx, double *y, int incy, 
		double c, double s)
{
  while (dim--){
    const double _x = *x;
    const double _y = *y;
    *x = c*_x + s*_y;
    *y = c*_y - s*_x;
    x += incx;
    y += incy;
  }
}


inline
void cblas_dswap(int dim, double *x, int incx, double *y, int incy)
{
  while (dim--){
    const double tmp = *x;
    *x = *y;
    *y = tmp;
    x += incx;
    y += incy;
  }
}


inline
void cblas_dcopy(int dim, const double *x, int incx, double *y, int incy)
{
  while(dim--){
    *y = *x;
    x += incx;
    y += incy;
  }
}


inline
void cblas_dscal(int dim, double alpha, double *x, int incx)
{
  while(dim--){
    *x *= alpha;
    x += incx;
  }
}


inline
void cblas_daxpy(int dim, double alpha, 
		 const double *x, int incx, double *y, int incy)
{
  while(dim--){
    *y += alpha * (*x);
    x += incx;
    y += incy;
  } 

  //while(dim--) y[incy * dim] += alpha * x[incx * dim];
}


inline 
double cblas_ddot(int dim, const double *x, int incx, 
		  const double *y, int incy)
{
  register double dot = 0.0;

  while(dim--){
    dot += (*x) * (*y);
    x += incx;
    y += incy;
  } 

  return dot;
}


inline
double cblas_dnrm2(int dim, const double *x, int incx)
{
  register double norm_sq = 0;

  while(dim--){
    norm_sq += (*x) * (*x);
    x += incx;
  } 

  return sqrt(norm_sq);
}


inline
int cblas_idamax(int dim, const double *x, int incx)
{
  double max_abs = fabs(*x);
  int argmax = 0;
  x += incx;

  for(int i=1; i<dim; i++){
    const double x_abs = fabs(*x);
    if (x_abs > max_abs){
      max_abs = x_abs;
      argmax = i;
    }
    x += incx;
  }

  return argmax;
}




// some BLAS Level 2 functions -- inline implementation


// A \in \R^{n \times m}
// x \in \R^m
// y \in \R^n
//
// computes y = beta y + alpha A x
inline
void cblas_dgemv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans,
		 int m, int n, double alpha, const double *A, int lda,
		 const double *x, int incx, 
		 double beta, double *y, int incy)
{
  // limited support!
  assert(order == CblasRowMajor);
  assert(trans == CblasNoTrans);
  
  while (n--){
    const double *_x = x;
    double tmp = 0.0;
    int j = m;

    while (j--){
      tmp += (*A) * (*_x);
      A++;
      _x += incx;
    }
    
    *y = beta*(*y) + alpha*tmp;
    y += incy;
  }
}







#endif
