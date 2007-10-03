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

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    double x = arg[0];
    double y = arg[1];

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

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }


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


template <int dim, class DomainField, class Field> 
class BenchMark_2 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
  const Field delta_;
  const Field sqrtDelta_;
  const Field x1_;
  const Field x2_;
public:  
  virtual ~BenchMark_2() {}
  BenchMark_2(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1e5)
    , sqrtDelta_( sqrt(delta_) )
    , x1_ ( 8. * atan(1.) )
    , x2_ ( x1_ / sqrtDelta_ ) 
  {
    factor_[0][0] = 1;
    factor_[0][1] = factor_[1][0] = 0.0;
    factor_[1][1] = delta_;
  }

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }


  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    return sin(x1_* x[0])*exp(-x2_* x[1]);
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = x1_ * cos(x1_* x[0])*exp(-x2_ * x[1]);
    grad[1] = -x2_ * exact(x);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    // we have only neumann boundary here 
    return false; 
  }
};



template <int dim, class DomainField, class Field> 
class BenchMark_3 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
  const Field delta_;
  const Field pi_;
  const Field cost_;
  const Field sint_;
public:  
  virtual ~BenchMark_3() {}
  BenchMark_3(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1e-3)
    , pi_ ( 4. * atan(1.) )
    , cost_ ( cos( 40. * pi_ / 180. ) )
    , sint_ ( sqrt(1. - cost_*cost_) ) 
  {
    DUNE_THROW(NotImplemented,"Not yet implemented");
    factor_[0][0] = cost_*cost_+delta_*sint_*sint_;
    factor_[1][0] = factor_[0][1] = cost_*sint_*(1.-delta_);
    factor_[1][1] = sint_*sint_+delta_*cost_*cost_;
  }

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    // u in general is unknown here
    // we only know u at the boundary 
    if( x[0] <  0.2 && x[1] <= 0.0 ) return 1;
    if( x[0] <= 0.0 && x[1] <  0.2 ) return 1;

    if( x[0] >  0.8 && x[1] >= 1.0 ) return 0;
    if( x[0] >= 1.0 && x[1] >  0.8 ) return 0;

    if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    if( x[0] <= 0.0 && x[1] >  0.3 ) return 0.5;

    if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    if( x[0] <  0.7 && x[1] >= 1.0 ) return 0.5;
    return 0.5;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 0.0;
    grad[1] = 0.0;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 

    if( x[0] > 0.2 && x[0] < 0.3 && x[1] <= 0.0) return false;
    if( x[0] > 0.7 && x[0] < 0.8 && x[1] >= 1.0) return false;

    if( x[0] >= 1.0 && x[1] > 0.7 && x[1] < 0.8) return false;
    if( x[0] <= 0.0 && x[1] > 0.2 && x[1] < 0.3) return false;
    // we have neumann boundary here 
    return true; 
  }
};


template <int dim, class DomainField, class Field> 
class BenchMark_4 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
  const Field delta_;
  const Field pi_;
  const Field cost_;
  const Field sint_;
public:  
  virtual ~BenchMark_4() {}
  BenchMark_4(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1)
    , pi_ ( 4. * atan(1.) )
    , cost_ ( cos( 40. * pi_ / 180. ) )
    , sint_ ( sqrt(1. - cost_*cost_) ) 
  {
    DUNE_THROW(NotImplemented,"Not yet implemented");
    factor_[0][0] = cost_*cost_+delta_*sint_*sint_;
    factor_[1][0] = factor_[0][1] = cost_*sint_*(1.-delta_);
    factor_[1][1] = sint_*sint_+delta_*cost_*cost_;
  }

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    // u in general is unknown here
    // we only know u at the boundary 
    if( x[0] <  0.2 && x[1] <= 0.0 ) return 1;
    if( x[0] <= 0.0 && x[1] <  0.2 ) return 1;

    if( x[0] >  0.8 && x[1] >= 1.0 ) return 0;
    if( x[0] >= 1.0 && x[1] >  0.8 ) return 0;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <= 0.0 && x[1] >  0.3 ) return 0.5;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <  0.7 && x[1] >= 1.0 ) return 0.5;
    return 0.5;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 0.0;
    grad[1] = 0.0;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 

    /*
    if( x[0] > 0.2 && x[0] < 0.3 && x[1] <= 0.0) return false;
    if( x[0] > 0.7 && x[0] < 0.8 && x[1] >= 1.0) return false;

    if( x[0] >= 1.0 && x[1] > 0.7 && x[1] < 0.8) return false;
    if( x[0] <= 0.0 && x[1] > 0.2 && x[1] < 0.3) return false;
*/
    // we have neumann boundary here 
    return true; 
  }
};


template <int dim, class DomainField, class Field> 
class BenchMark_5 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field delta_;
  const Field pi;
public:  
  virtual ~BenchMark_5() {}
  BenchMark_5(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1e-3)
    , pi ( 4. * atan(1.) )
  {
  }

  virtual void factor(const DomainField arg[dim], Field k[dim][dim] ) const
  {
    double x = arg[0];
    double y = arg[1];
    double rt = x*x+y*y;
    k[0][0] = (y*y+delta_*x*x)/rt;
    k[1][1] = (x*x+delta_*y*y)/rt;
    k[1][0] = k[0][1] = -(1-delta_)*x*y/rt;
  }

  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    Field k[dim][dim];
    factor(arg,k);
    double x = arg[0];
    double y = arg[1];
    double rt = x*x+y*y;

    double ux = pi * cos(pi*x)*sin(pi*y);
    double uy = pi * cos(pi*y)*sin(pi*x);

    double f0 = sin(pi*x)*sin(pi*y)*pi*pi*(1+delta_)*(x*x+y*y) 
              + cos(pi*x)*sin(pi*y)*pi*(1.-3.*delta_)*x
              + cos(pi*y)*sin(pi*x)*pi*(1.-3.*delta_)*y
              + cos(pi*y)*cos(pi*x)*2.*pi*pi*(1.-delta_)*x*y;
    double kxx = k[0][0];
    double kyy = k[1][1];
    double kxy = k[0][1];
    return (f0+2.*(x*(kxx*ux+kxy*uy)+y*(kxy*ux+kyy*uy)))/rt;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    return sin(pi*x[0])*sin(pi*x[1]);
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = pi * cos(pi*x[0])*sin(pi*x[1]);
    grad[1] = pi * cos(pi*x[1])*sin(pi*x[0]);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


template <int dim, class DomainField, class Field> 
class BenchMark_6 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
  const Field delta_;
  const Field pi_;
  const Field cost_;
  const Field sint_;
public:  
  virtual ~BenchMark_6() {}
  BenchMark_6(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1)
    , pi_ ( 4. * atan(1.) )
    , cost_ ( cos( 40. * pi_ / 180. ) )
    , sint_ ( sqrt(1. - cost_*cost_) ) 
  {
    DUNE_THROW(NotImplemented,"Not yet implemented");
    factor_[0][0] = cost_*cost_+delta_*sint_*sint_;
    factor_[1][0] = factor_[0][1] = cost_*sint_*(1.-delta_);
    factor_[1][1] = sint_*sint_+delta_*cost_*cost_;
  }

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }


  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    // u in general is unknown here
    // we only know u at the boundary 
    if( x[0] <  0.2 && x[1] <= 0.0 ) return 1;
    if( x[0] <= 0.0 && x[1] <  0.2 ) return 1;

    if( x[0] >  0.8 && x[1] >= 1.0 ) return 0;
    if( x[0] >= 1.0 && x[1] >  0.8 ) return 0;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <= 0.0 && x[1] >  0.3 ) return 0.5;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <  0.7 && x[1] >= 1.0 ) return 0.5;
    return 0.5;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 0.0;
    grad[1] = 0.0;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 

    /*
    if( x[0] > 0.2 && x[0] < 0.3 && x[1] <= 0.0) return false;
    if( x[0] > 0.7 && x[0] < 0.8 && x[1] >= 1.0) return false;

    if( x[0] >= 1.0 && x[1] > 0.7 && x[1] < 0.8) return false;
    if( x[0] <= 0.0 && x[1] > 0.2 && x[1] < 0.3) return false;
*/
    // we have neumann boundary here 
    return true; 
  }
};


template <int dim, class DomainField, class Field> 
class BenchMark_7 : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  Field factor_[dim][dim];
  const Field delta_;
  const Field pi_;
  const Field cost_;
  const Field sint_;
public:  
  virtual ~BenchMark_7() {}
  BenchMark_7(Field globalShift, Field factor)
    : globalShift_(0.0)
    , delta_(1)
    , pi_ ( 4. * atan(1.) )
    , cost_ ( cos( 40. * pi_ / 180. ) )
    , sint_ ( sqrt(1. - cost_*cost_) ) 
  {
    DUNE_THROW(NotImplemented,"Not yet implemented");
    factor_[0][0] = cost_*cost_+delta_*sint_*sint_;
    factor_[1][0] = factor_[0][1] = cost_*sint_*(1.-delta_);
    factor_[1][1] = sint_*sint_+delta_*cost_*cost_;
  }

  virtual void factor(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i)
    {
      k[i][i] = factor_[i][i];
      for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
      for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
    }
  }


  virtual Field rhs  (const DomainField arg[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    // u in general is unknown here
    // we only know u at the boundary 
    if( x[0] <  0.2 && x[1] <= 0.0 ) return 1;
    if( x[0] <= 0.0 && x[1] <  0.2 ) return 1;

    if( x[0] >  0.8 && x[1] >= 1.0 ) return 0;
    if( x[0] >= 1.0 && x[1] >  0.8 ) return 0;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <= 0.0 && x[1] >  0.3 ) return 0.5;

    //if( x[0] >  0.3 && x[1] <= 0.0 ) return 0.5;
    //if( x[0] <  0.7 && x[1] >= 1.0 ) return 0.5;
    return 0.5;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 0.0;
    grad[1] = 0.0;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 

    /*
    if( x[0] > 0.2 && x[0] < 0.3 && x[1] <= 0.0) return false;
    if( x[0] > 0.7 && x[0] < 0.8 && x[1] >= 1.0) return false;

    if( x[0] >= 1.0 && x[1] > 0.7 && x[1] < 0.8) return false;
    if( x[0] <= 0.0 && x[1] > 0.2 && x[1] < 0.3) return false;
*/
    // we have neumann boundary here 
    return true; 
  }
};
#endif
