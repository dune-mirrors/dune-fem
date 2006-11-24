#ifndef SUBBLAS_HPP
#define SUBBLAS_HPP

#include <cmath>
#include <cstring>
#include <cassert>


// a subset of BLAS subroutines


// Enumerated and derived types
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};


// BLAS level 1 functions
void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drot(const int dim, double *x, const int incx, 
		double *y, const int incy, const double c, const double s);
void cblas_dswap(const int dim, double *x, const int incx, 
		 double *y, const int incy);
void cblas_dcopy(const int dim, const double *x, const int incx, 
		 double *y, const int incy);
void cblas_dscal(const int dim, const double alpha, double *x, const int incx);
void cblas_daxpy(const int dim, const double alpha, 
		 const double *x, const int incx, double *y, const int incy);
double cblas_ddot(const int dim, const double *x, const int incx, 
		  const double *y, const int incy);
double cblas_dnrm2(const int dim, const double *x, const int incx);
int cblas_idamax(const int dim, const double *x, const int incx);


// BLAS level 2 functions
void cblas_dgemv(const enum CBLAS_ORDER order, 
		 const enum CBLAS_TRANSPOSE trans,
		 const int m, const int n, 
		 const double alpha, const double *A, const int lda,
		 const double *x, const int incx, 
		 const double beta, double *y, const int incy);


// BLAS level 3 functions
void cblas_dgemm(const enum CBLAS_ORDER order, 
		 const enum CBLAS_TRANSPOSE transA, 
		 const enum CBLAS_TRANSPOSE transB, 
		 const int m, const int n, const int k, 
		 const double alpha, const double *A, const int lda, 
		 const double *B, const int ldb,
                 const double beta, double *C, const int ldc);



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
void cblas_drot(const int dim, double *x, const int incx, 
		double *y, const int incy, 
		double c, double s)
{
  int i = dim;

  while (i--){
    const double _x = *x;
    const double _y = *y;
    *x = c*_x + s*_y;
    *y = c*_y - s*_x;
    x += incx;
    y += incy;
  }
}


inline
void cblas_dswap(const int dim, double *x, const int incx, 
		 double *y, const int incy)
{
  int i = dim;

  while (i--){
    const double tmp = *x;
    *x = *y;
    *y = tmp;
    x += incx;
    y += incy;
  }
}


inline
void cblas_dcopy(const int dim, const double *x, const int incx, 
		 double *y, const int incy)
{
  if (incx == 1 && incy == 1){ // use memcpy
    memcpy(y, x, dim*sizeof(double));
  }
  else{
    int i = dim;
    while(i--){
      *y = *x;
      x += incx;
      y += incy;
    }
  }
}


inline
void cblas_dscal(const int dim, const double alpha, double *x, const int incx)
{
  int i = dim;
  while(i--){
    *x *= alpha;
    x += incx;
  }
}


inline
void cblas_daxpy(const int dim, const double alpha, 
		 const double *x, const int incx, double *y, const int incy)
{
  int i = dim;
  while(i--){
    *y += alpha * (*x);
    x += incx;
    y += incy;
  } 
}


inline 
double cblas_ddot(const int dim, const double *x, const int incx, 
		  const double *y, const int incy)
{
  register double dot = 0.0;
  int i=dim;

  while(i--){
    dot += (*x) * (*y);
    x += incx;
    y += incy;
  } 

  return dot;
}


inline
double cblas_dnrm2(const int dim, const double *x, const int incx)
{
  register double norm_sq = 0;
  int i = dim;

  while(i--){
    norm_sq += (*x) * (*x);
    x += incx;
  } 

  return sqrt(norm_sq);
}


inline
int cblas_idamax(const int dim, const double *x, const int incx)
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


// op(A) \in \R^{m \times n}
// x \in \R^n
// y \in \R^m
//
// computes y = beta y + alpha op(A) x
inline
void cblas_dgemv(const enum CBLAS_ORDER order, 
		 const enum CBLAS_TRANSPOSE trans,
		 const int m, const int n, 
		 const double alpha, const double *A, const int lda,
		 const double *x, const int incx, 
		 const double beta, double *y, const int incy)
{
  // limited support!
  assert(order == CblasRowMajor);

  if (trans == CblasNoTrans){
    assert(lda >= n); // consistent lda
  
    for(int i=0; i<m; i++){
      double sum = 0.0;
      for(int j=0; j<n; j++) sum += A[lda*i + j] * x[incx*j];
      y[incy*i] = beta * y[incy*i] + alpha * sum;
    }     
  }
  else{
    assert(trans == CblasTrans);
    assert(lda >= m);  // consistent lda  

    for(int i=0; i<m; i++){
      double sum = 0.0;
      for(int j=0; j<n; j++) sum += A[lda*j + i] * x[incx*j];
      y[incy*i] = beta * y[incy*i] + alpha * sum;
    }     
  }
}




// some BLAS Level 3 functions -- inline implementation


// A \in \R^{m \times k}, lda>=k
// B \in \R^{k \times n}, ldb>=n
// C \in \R^{m \times n}, ldc>=n
//
// computes C = alpha A*B + beta C
inline
void cblas_dgemm(const enum CBLAS_ORDER order, 
		 const enum CBLAS_TRANSPOSE transA, 
		 const enum CBLAS_TRANSPOSE transB, 
		 const int m, const int n, const int k, 
		 const double alpha, const double *A, const int lda, 
		 const double *B, const int ldb,
                 const double beta, double *C, const int ldc)
{
  // limited support!
  assert(order == CblasRowMajor);
  assert(transB == CblasNoTrans);

  // new
  if (transA == CblasNoTrans){
    // consistent lda, ldb, ldc
    assert(lda >= k && ldb >= n && ldc >= n);

    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
	double sum = 0.0;
 	for(int l=0; l<k; l++) sum += A[i*lda + l] * B[l*ldb + j];
 	C[i*ldc + j] = alpha*sum + beta*C[i*ldc + j];
      }
    }
  }
  else{
    assert (transA == CblasTrans);

    // consistent lda, ldb, ldc
    assert(lda >= m && ldb >= n && ldc >= n);

    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
	double sum = 0.0;
 	for(int l=0; l<k; l++) sum += A[l*lda + i] * B[l*ldb + j];
 	C[i*ldc + j] = alpha*sum + beta*C[i*ldc + j];
       }
    }
  }
}






#endif
