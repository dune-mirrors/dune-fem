#ifndef FVCA5_BENCHMARK_CC
#define FVCA5_BENCHMARK_CC

#include <cmath>
#include "problem.cc"

template <int dim, class DomainField, class Field> 
class BenchMark_1 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
public:  
  virtual ~BenchMark_1() {}
  BenchMark_1(Field globalShift, Field factor)
    : globalShift_(0.0)
  {
    for(int i=0; i<dim; ++i) 
    {
      for(int j=0; j<dim; ++j) 
      {
        if( i == j ) factor_[i][j] = 1.5;
        else factor_[i][j] = 0.5;
      }
    }
  }

  virtual Field factor (int i, int j) const { return factor_[i][j]; }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    double x = arg[0];
    double y = arg[1];

    //u  = x*(1-x)*y*(1-y)*16.

    //double ux = (-2.*x+1)*y*(1-y)*16.
    //double uy = (-2.*y+1)*x*(1-x)*16.

    double uxx = -2.*y*(1-y)*16.;
    double uxy = (-2.*x+1)*(-2.*y+1)*16.;
    double uyy = -2.*x*(1-x)*16.;

    return -( factor_[0][0] * uxx + factor_[0][1] * uxy + factor_[1][0] * uxy + factor_[1][1] * uyy);
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = 16.;
    for(int i=0; i<dim; ++i) 
      val *= (x[i] - (x[i]*x[i]));
    return val;
  }
  
  virtual void gradExact(const DomainField arg[dim], Field grad[dim] ) const 
  {
    double x = arg[0];
    double y = arg[1];

    grad[0] = (-2.*x+1)*y*(1-y)*16.;
    grad[1] = (-2.*y+1)*x*(1-x)*16.;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


template <int dim, class DomainField, class Field> 
class BenchMark_1_2 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
public:  
  virtual ~BenchMark_1_2() {}
  BenchMark_1_2(Field globalShift, Field factor)
    : globalShift_(0.0)
  {
    for(int i=0; i<dim; ++i) 
    {
      for(int j=0; j<dim; ++j) 
      {
        if( i == j ) factor_[i][j] = 1.5;
        else factor_[i][j] = 0.5;
      }
    }
  }

  virtual Field factor (int i, int j) const { return factor_[i][j]; }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    double x1 = 1.-arg[0];
    double y1 = 1.-arg[1];

    double uxx = -y1*y1*sin(x1*y1) + 6.*x1* y1*y1;
    double uxy = -x1*y1*sin(x1*y1) + cos(x1*y1) + 6.*y1* x1*x1;
    double uyy = -x1*x1*sin(x1*y1) + 2.* x1*x1*x1;
    return -(factor_[0][0] * uxx + factor_[0][1]*uxy + factor_[1][0]* uxy + factor_[1][1]*uyy);
  }

  virtual Field exact(const DomainField arg[dim]) const
  {
    double x1 = 1.-arg[0];
    double y1 = 1.-arg[1];
    return sin(x1*y1) + x1*x1*x1 * y1*y1;
  }
  
  virtual void gradExact(const DomainField arg[dim], Field grad[dim] ) const 
  {
    double x1 = 1.-arg[0];
    double y1 = 1.-arg[1];

    grad[0] = -y1*cos(x1*y1) - 3.* (x1*y1)* (x1*y1);
    grad[1] = -x1*cos(x1*y1) - 2.*y1* x1*x1*x1;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};



#endif
