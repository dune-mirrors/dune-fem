#ifndef ELLIPTPROBLEM_CC
#define ELLIPTPROBLEM_CC

#include <cmath>

template<int dim> 
class DataFunctionIF
{
protected:  
  DataFunctionIF() {}
public:  
  virtual ~DataFunctionIF() {}
  // factor 
  virtual double factor() const = 0;
  // right hand side 
  virtual double rhs  (const double x[dim]) const = 0;
  // exact solution 
  virtual double exact(const double x[dim]) const = 0;
  // exact gradient 
  virtual void gradExact(const double x[dim], double grad[dim] ) const = 0;

  // boundary data 
  virtual bool boundaryDataFunction(const double x[dim], double & val) const = 0;

  // neumann boundary 
  virtual void neumann(const double x[dim], double grad[dim]) const 
  {
    gradExact(x,grad);
    for(int i=0; i<dim; ++i) grad[i] *= factor();
  }
};


template <int dim> 
class SinSin : public DataFunctionIF<dim>
{
  const double globalShift_;
  const double factor_;
public:  
  virtual ~SinSin() {}
  SinSin(double globalShift, double factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    //assert(dim == 2);
  }

  virtual double factor () const { return factor_; }

  virtual double rhs  (const double x[dim]) const 
  {
    double sin_x = sin(2.0*M_PI*x[0]);
    double sin_y = sin(2.0*M_PI*x[1]);

    double val = 8.0 * M_PI * M_PI * sin_x * sin_y ;
    val *= factor_;
    return val;
  }

  virtual double exact(const double x[dim]) const
  {
    double val = sin(2.0*M_PI*x[0]) * sin(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const double x[dim], double grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*cos(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
  }
  
  virtual bool boundaryDataFunction(const double x[dim], double & val) const
  {
    val = exact( x ); 
    if(x[0] <= 0.0) return false;
    return true; 
  }
};

template <int dim> 
class CosCos : public DataFunctionIF<dim>
{
  const double globalShift_;
  const double factor_;
public:  
  virtual ~CosCos() {}
  CosCos(double globalShift, double factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    //assert(dim == 2);
  }

  virtual double factor () const { return factor_; }

  virtual double rhs  (const double x[dim]) const 
  {
    double cos_x = cos(2.0*M_PI*x[0]);
    double cos_y = cos(2.0*M_PI*x[1]);
       
    double val = 8.0 * M_PI*M_PI* cos_x * cos_y ;
    val *= factor_;
    return val;
  }

  virtual double exact(const double x[dim]) const
  {
    double val = cos(2.0*M_PI*x[0]) * cos(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const double x[dim], double grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*-sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*cos(2.0*M_PI*x[0])*-sin(2.0*M_PI*x[1]);
  }
  
  virtual bool boundaryDataFunction(const double x[dim], double & val) const
  {
    val = exact( x ); 
    if(x[0] <= 0.0) return false;
    return true; 
  }
};

//! problem from Castillo paper 
template <int dim> 
class CastilloProblem : public DataFunctionIF<dim>
{
  const double globalShift_;
  const double factor_;
public:  
  virtual ~CastilloProblem() {}
  CastilloProblem(double globalShift, double factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    assert(dim == 2);
  }

  virtual double factor () const { return factor_; }

  virtual double rhs  (const double x[dim]) const 
  {
    double ret = 0.0;
    double tmp =-23+7*SQR(x[1])-24*x[0]+24*x[0]*SQR(x[1])+7*SQR(x[0])+9*SQR(x[0])*SQR(x[1])-24
                *x[1]+24*x[1]*SQR(x[0]);

    ret=-0.5*exp(0.75*(x[0]+x[1]));
    ret *= tmp;
    ret *= factor_;
    return ret;
  }

  virtual double exact(const double x[dim]) const
  {
    double val = 4.0 *(1.-SQR(x[0]))*(1.-SQR(x[1]))*exp(0.75*(x[0]+x[1]));
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const double x[dim], double grad[dim] ) const 
  {
    // not implemented yet 
    grad[0] = grad[1] = 0.0; 
  }
  
  virtual bool boundaryDataFunction(const double x[dim], double & val) const
  {
    val = exact( x ); 
    return true; 
  }
};

//! problem Einstringende Ecke 
template <int dim> 
class InSpringingCorner : public DataFunctionIF<dim>
{
  const double globalShift_;
  const double factor_;
public:  
  virtual ~InSpringingCorner() {}
  InSpringingCorner(double globalShift, double factor)
    : globalShift_(0.0)
    , factor_(1.0) 
  {
    assert(dim == 2);
  }

  virtual double factor () const { return factor_; }

  virtual double rhs  (const double x[dim]) const 
  {
    return 0.0;
  }

  virtual double exact(const double x[dim]) const
  {
    double ret = 0.0;
    const double Pi= M_PI; 
    double tmp2=0.0;

    double s = 0.0 ;
    for(int i=0; i<dim; ++i) s += x[i]*x[i];
    double r = sqrt(s); 
   
    double phi=0.0;
    double fac=2.0;
    fac/=3;
   
    ret=pow(r,fac);
    if(x[1]>=0.0)
    {
      tmp2=x[0];
      tmp2/=r;
      phi=acos(tmp2);
      phi*=fac;
      ret*=sin(phi);
    }
    else
    {
      tmp2=-x[0];
      tmp2/=r;
      phi=acos(tmp2);
      phi+=Pi;
      phi*=fac;
      ret*=sin(phi);
    }
    return ret;
  }
  
  virtual void gradExact(const double x[dim], double grad[dim] ) const 
  {
    // not implemented yet 
    grad[0] = grad[1] = 0.0; 
  }
  
  virtual bool boundaryDataFunction(const double x[dim], double & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


//! problem Riviere-Bastian 
template <int dim> 
class RiviereProblem : public DataFunctionIF<dim>
{
  const double globalShift_;
  const double factor_;
public:  
  virtual ~RiviereProblem() {}
  RiviereProblem(double globalShift, double factor)
    : globalShift_(0.0)
    , factor_(1.0) 
  {
    assert(dim == 2);
  }

  virtual double factor () const { return factor_; }

  virtual double rhs  (const double x[dim]) const 
  {
    double val = exact(x);
    double x_part = val * SQR(-2.0 * (x[0] - 0.5)) - 2.0 * val; 
    double y_part = val * SQR(-2.0 * (x[1] - 0.5)) - 2.0 * val; 
    return -(x_part + y_part);
  }

  virtual double exact(const double x[dim]) const
  {
    double power = -( SQR( x[0] - 0.5 ) + SQR( x[1] - 0.5 ) );
    return pow( M_E , power );
  }
  
  virtual void gradExact(const double x[dim], double grad[dim] ) const 
  {
    double val = exact( x );
    grad[0] = val * ( -2.0*x[0] + 1.0 );
    grad[1] = val * ( -2.0*x[1] + 1.0 );
  }
  
  virtual bool boundaryDataFunction(const double x[dim], double & val) const
  {
    val = exact( x ); 
    return true; 
  }
};
#endif
