#include "matrix.hpp"
#include "linear_solver.hpp"


//using namespace pardg;


double pardg::Matrix::output_epsilon = 1.0e-7;



pardg::Matrix& pardg::Matrix::inverse()
{
  assert(n == m);
  transpose();
  QRSolver linear_solver;
  linear_solver.prepare(m, data);

  double *data_inv = new double[n*n];
  double *_data_inv = data_inv;
  for(int i=0; i<n; i++){
    dset(n, 0.0, _data_inv, 1);
    _data_inv[i] = 1.0;
    linear_solver.solve(_data_inv);
    _data_inv += n;
  }

  cblas_dcopy(n*n, data_inv, 1, data, 1);   
  delete[] data_inv;
  return *this;
}



pardg::Matrix& pardg::Matrix::transpose()
{
  assert(n == m);
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      double tmp = data[i*n + j];
      data[i*n + j] = data[j*n + i];
      data[j*n + i] = tmp;      
    }
  }
  return *this;
}



pardg::Matrix& pardg::Matrix::identity()
{
  assert(n == m);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      data[i*n + j] = (i==j)? 1.0: 0.0;
    }
  }
  return *this;
}



pardg::Matrix pardg::operator+(const pardg::Matrix &A, const pardg::Matrix &B)
{
  assert(A.n == B.n && A.m == B.m);
  Matrix C(A.n, A.m);
  dwaxpby(A.n*A.m, 1.0, A, 1, 1.0, B, 1, C, 1);
  return C;
}



pardg::Matrix pardg::operator*(const pardg::Matrix &A, const pardg::Matrix &B)
{
  assert(A.m == B.n);
  const int n = A.n;
  const int m = B.m;
  Matrix C(n, m);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double sum = 0.0;
      for(int k=0; k<A.m; k++) sum += A(i,k) * B(k,j);
      C(i,j) = sum;
    }
  }
    
  return C;
}



pardg::Matrix pardg::operator*(double lambda, const pardg::Matrix &A)
{
  Matrix B(A.n, A.m);
  cblas_daxpy(A.n*A.m, lambda, A, 1, B, 1);
  return B;
}




std::ostream& pardg::operator<<(std::ostream& os, const pardg::Matrix& A)
{
  for(int i=0; i<A.n; i++){
    for(int j=0; j<A.m; j++){
      os.width(12);
      if ( fabs(A(i,j)) >= Matrix::output_epsilon ) os << A(i,j) << "  ";
      else  os << 0.0 << "  ";
    }
    os << std::endl;
  }
  return os;
}

