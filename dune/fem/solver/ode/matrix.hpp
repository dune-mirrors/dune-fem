#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cassert>
#include "function.hpp"
#include "blas.hpp"


namespace pardg
{


class Matrix : public Function
{
public:
  Matrix(int n, int m);
  Matrix(int n);
  Matrix(const Matrix &A); // copy constructor
  Matrix(int n, int m, const double *a);
  Matrix(int n, const double *a);   
  virtual ~Matrix();
  
  // element access
  double& operator()(int i, int j);
  double operator()(int i, int j) const;

  // conversion operators
  operator double*();
  operator const double *() const;

  // assignment operators
  Matrix& operator=(const Matrix& A);
  Matrix& operator=(const double *a);
  Matrix& operator+=(const Matrix& A);
  Matrix& operator-=(const Matrix& A);
  Matrix& operator*=(double alpha);
  Matrix& operator=(double alpha);

  // other operations for (n x n) matrices
  Matrix& inverse();
  Matrix& transpose();
  Matrix& identity();

  // from Function
  virtual void operator()(const double *u, double *f, int i=0);
  virtual int dim_of_argument(int i=0) const;
  virtual int dim_of_value(int i=0) const;

  // not efficient but sometimes useful operators
  friend Matrix operator+(const Matrix &A, const Matrix &B);
  friend Matrix operator*(const Matrix &A, const Matrix &B);
  friend Matrix operator*(double lambda, const Matrix &A);
  friend std::ostream& operator<<(std::ostream& os, const Matrix& A);

  static double output_epsilon;
  
private:
  const int n, m;
  double *data;
};


} // namespace pardg




inline 
pardg::Matrix::Matrix(int n, int m) : n(n), m(m), data(new double[n*m])
{
  assert(data);
  dset(n*m, 0.0, data, 1);
}


inline 
pardg::Matrix::Matrix(int n) : n(n), m(n), data(new double[n*n])
{
  assert(data);
  dset(n*m, 0.0, data, 1);
}


inline 
pardg::Matrix::Matrix(const Matrix &A) : n(A.n), m(A.m), data(new double[A.n*A.m])
{
  assert(data);
  cblas_dcopy(n*m, A.data, 1, data, 1);
}


inline 
pardg::Matrix::Matrix(int n, int m, const double *a) : 
  n(n), m(m), data(new double[n*m])
{
  assert(data);
  cblas_dcopy(n*m, a, 1, data, 1);
}


inline 
pardg::Matrix::Matrix(int n, const double *a) : 
  n(n), m(n), data(new double[n*n])
{
  assert(data);
  cblas_dcopy(n*m, a, 1, data, 1);
}


inline 
pardg::Matrix::~Matrix()
{
  delete[] data;
}


inline
double& pardg::Matrix::operator()(int i, int j)
{
  assert(i>=0 && i<n && j>=0 && j<m);
  return data[i*m + j];
}


inline
double pardg::Matrix::operator()(int i, int j) const
{
  assert(i>=0 && i<n && j>=0 && j<m);
  return data[i*m + j];
}


inline
pardg::Matrix::operator double*()
{
  return data;
}


inline
pardg::Matrix::operator const double *() const
{
  return data;
}


inline
pardg::Matrix& pardg::Matrix::operator=(const Matrix& A)
{
  assert(n==A.n && m==A.m);
  cblas_dcopy(n*m, A.data, 1, data, 1);
  return *this;
}


inline
pardg::Matrix& pardg::Matrix::operator=(const double *a)
{
  cblas_dcopy(n*m, a, 1, data, 1);
  return *this;
}


inline
pardg::Matrix& pardg::Matrix::operator+=(const Matrix& A)
{
  assert(n==A.n && m==A.m);
  cblas_daxpy(n*m, 1.0, A.data, 1, data, 1);
  return *this;
}


inline
pardg::Matrix& pardg::Matrix::operator-=(const Matrix& A)
{
  assert(n==A.n && m==A.m);
  cblas_daxpy(n*m, -1.0, A.data, 1, data, 1);
  return *this;
}


inline
pardg::Matrix& pardg::Matrix::operator*=(double alpha)
{
  cblas_dscal(n*m, alpha, data, 1);
  return *this;
}


inline
pardg::Matrix& pardg::Matrix::operator=(double alpha)
{
  dset(n*m, alpha, data, 1);
  return *this;
}


inline
void pardg::Matrix::operator()(const double *u, double *f, int i)
{
  switch (i){
  case 0: // function
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
		m, n, 1.0, data, m, u, 1, 0.0, f, 1);  
    break;

  case 1: // first derivative
    cblas_dcopy(n*m, data, 1, f, 1);
    break;

  default:
    assert(0);
  } 
}


inline 
int pardg::Matrix::dim_of_argument(int i) const
{
  switch (i){
  case 0: return m;
  case 1: return n*m;
  default: assert(0); return -1;
  }
}


inline
int pardg::Matrix::dim_of_value(int i) const
{
  switch (i){
  case 0: return n;
  case 1: return n*m;
  default: assert(0); return -1;
  }
}



#endif
