#ifndef FEM_BLASVECTOR_HH
#define FEM_BLASVECTOR_HH

#include <iostream>
#include <cassert>
#include <cmath>
#include <dune/fem/solver/pardg.hh>
#include <dune/fem/storage/vector.hh>
namespace Dune {

class BlasVector : 
public VectorDefault<double,BlasVector> {
  bool owner_;
public:
  typedef double FieldType;
  BlasVector(unsigned int n);
  BlasVector(unsigned int n, const double *a);
  BlasVector(const BlasVector &v); // copy constructor
  ~BlasVector();
  
  // element access
  double& operator[](unsigned int i);
  double operator[](unsigned int i) const;

  // assignment operators
  BlasVector& operator=(const BlasVector& v);
  BlasVector& operator+=(const BlasVector& v);
  BlasVector& operator-=(const BlasVector& v);
  BlasVector& operator*=(double alpha);
  BlasVector& operator=(double alpha);

  unsigned int size() const;

  double* leakPointer() {
    return data;
  }
  const double* leakPointer() const {
    return data;
  }
  void reserve(unsigned int newSize) {
    assert(owner_);
    if (newSize<n) return;
    unsigned int useSize = newSize;
    double* newData = new double[useSize];
    if (!data) data = newData;
    else {
      double* oldData;
      pardg::cblas_dcopy(n,data,1,newData,1);
      data = newData;
      delete [] oldData;
    }
    totalSize = useSize;

  }
  void resize(unsigned int newSize) {
    assert(owner_);
    if (newSize>totalSize) 
      reserve(newSize);
    else 
      n = newSize;
  }
private:
  unsigned int n,totalSize;
  double *data;
};


// class BlasVector inline implementation

inline 
BlasVector::BlasVector(unsigned int pn) : owner_(true),
n(pn), totalSize(pn), data(new double[pn])
{
  assert(data);
  pardg::dset(n, 0.0, data, 1);
}


inline 
BlasVector::BlasVector(unsigned int pn, const double *a) : owner_(false),
n(pn), totalSize(pn), data(const_cast<double*>(a))
{
  assert(data);
  // pardg::cblas_dcopy(n, a, 1, data, 1);
}


inline 
BlasVector::BlasVector(const BlasVector &v) : owner_(true),
n(v.n), totalSize(v.n), data(new double[v.n])
{
  assert(data);
  pardg::cblas_dcopy(n, v.data, 1, data, 1);
}


inline 
BlasVector::~BlasVector()
{
  if (owner_)
    delete[] data;
}

inline
double& BlasVector::operator[](unsigned int i)
{
  assert(i<n);
  return data[i];
}


inline
double BlasVector::operator[](unsigned int i) const
{
  assert(i>=0 && i<n);
  return data[i];
}


inline
BlasVector& BlasVector::operator=(const BlasVector& v)
{
  assert(n==v.n);
  pardg::cblas_dcopy(n, v.data, 1, data, 1);
  return *this;
}


inline
BlasVector& BlasVector::operator+=(const BlasVector& v)
{
  assert(n==v.n);
  pardg::cblas_daxpy(n, 1.0, v.data, 1, data, 1);
  return *this;
}


inline
BlasVector& BlasVector::operator-=(const BlasVector& v)
{
  assert(n==v.n);
  pardg::cblas_daxpy(n, -1.0, v.data, 1, data, 1);
  return *this;
}


inline
BlasVector& BlasVector::operator*=(double alpha)
{
  pardg::cblas_dscal(n, alpha, data, 1);
  return *this;
}


inline
BlasVector& BlasVector::operator=(double alpha)
{
  pardg::dset(n, alpha, data, 1);
  return *this;
}


inline
unsigned int BlasVector::size() const
{
  return n;
}

}

#endif
