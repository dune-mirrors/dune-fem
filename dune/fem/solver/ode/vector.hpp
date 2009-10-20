#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <iostream>
#include <cassert>
#include <cmath>
#include "blas.hpp"


class Vector
{
public:
  Vector(int n);
  Vector(int n, const double *a);
  Vector(const Vector &v); // copy constructor
  ~Vector();
  
  // element access
  double& operator()(int i);
  double operator()(int i) const;
  double& operator[](int i);
  double operator[](int i) const;

  // conversion operators
  operator double*();
  operator const double *() const;

  // assignment operators
  Vector& operator=(const Vector& v);
  Vector& operator+=(const Vector& v);
  Vector& operator-=(const Vector& v);
  Vector& operator*=(double alpha);
  Vector& operator=(double alpha);

  int size() const;

  friend std::ostream& operator<<(std::ostream& os, const Vector& v);

  static double output_epsilon;

private:
  const int n;
  double *data;
};


// class Vector inline implementation

inline 
Vector::Vector(int n) : n(n), data(new double[n])
{
  assert(data);
  dset(n, 0.0, data, 1);
}


inline 
Vector::Vector(int n, const double *a) : n(n), data(new double[n])
{
  assert(data);
  cblas_dcopy(n, a, 1, data, 1);
}


inline 
Vector::Vector(const Vector &v) : n(v.n), data(new double[v.n])
{
  assert(data);
  cblas_dcopy(n, v.data, 1, data, 1);
}


inline 
Vector::~Vector()
{
  delete[] data;
}


inline
double& Vector::operator()(int i)
{
  assert(i>=0 && i<n);
  return data[i];
}


inline
double Vector::operator()(int i) const
{
  assert(i>=0 && i<n);
  return data[i];
}


inline
double& Vector::operator[](int i)
{
  assert(i>=0 && i<n);
  return data[i];
}


inline
double Vector::operator[](int i) const
{
  assert(i>=0 && i<n);
  return data[i];
}


inline
Vector::operator double*()
{
  return data;
}


inline
Vector::operator const double *() const
{
  return data;
}


inline
Vector& Vector::operator=(const Vector& v)
{
  assert(n==v.n);
  cblas_dcopy(n, v.data, 1, data, 1);
  return *this;
}


inline
Vector& Vector::operator+=(const Vector& v)
{
  assert(n==v.n);
  cblas_daxpy(n, 1.0, v.data, 1, data, 1);
  return *this;
}


inline
Vector& Vector::operator-=(const Vector& v)
{
  assert(n==v.n);
  cblas_daxpy(n, -1.0, v.data, 1, data, 1);
  return *this;
}


inline
Vector& Vector::operator*=(double alpha)
{
  cblas_dscal(n, alpha, data, 1);
  return *this;
}


inline
Vector& Vector::operator=(double alpha)
{
  dset(n, alpha, data, 1);
  return *this;
}


inline
int Vector::size() const
{
  return n;
}


inline
std::ostream& operator<<(std::ostream& os, const Vector& v)
{
  for(int i=0; i<v.n; i++){
    os.width(10);
    if (fabs(v(i)) >= Vector::output_epsilon) os << v(i) << std::endl;
    else os << 0.0 << std::endl;
  }
  return os;
}



#endif
