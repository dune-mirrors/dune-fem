#ifndef ELLIPTPROBLEM_CC
#define ELLIPTPROBLEM_CC

#include <math.h>

#define PROBLEM 1

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
    assert(dim == 2);
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
    assert(dim == 2);
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

#if 0
// shift solutions to avoid dirichlet 0 bnd 
static const double globalShift = 2.0;

double exactFactor() 
{
  return 1.0/(5.24*1e5);
  //return 12.0;
  //return 1.;
  //return (5.24*1e5);
}

double frhs(const double x[dim]) 
{
  double ret = 0.0;
  for(int i=0; i<dim; i++)
  {
    double tmp = 2.0;
    for(int j=1; j<dim; j++)
    {
      int idx = (i+j) % dim;
      tmp *= (x[idx] - SQR(x[idx]));
    } 
    ret += tmp;
  } 
  return ret;
} 

double coscosHalf(const double x[dim]) 
{
  double cos_x = cos(0.5*M_PI*x[0]);
  double cos_y = cos(0.5*M_PI*x[1]);
       
  double val = -0.5 * M_PI*M_PI* cos_x * cos_y ;
  return -val;
}

double expProblemRhs( const double x[dim]) 
{
  double ret = 0.0;
  double tmp =-23+7*SQR(x[1])-24*x[0]+24*x[0]*SQR(x[1])+7*SQR(x[0])+9*SQR(x[0])*SQR(x[1])-24
              *x[1]+24*x[1]*SQR(x[0]);

  ret=-0.5*exp(0.75*(x[0]+x[1]));
  ret*=tmp;
  return ret;
}

double rhsFunction(const double x[dim])
{
  double val = 0.0;
#if PROBLEM == 3
  val = M_PI * M_PI * cos( M_PI * x[0] );
#endif

#if PROBLEM == 4 
  val = frhs( x );
#endif
  
#if PROBLEM == 5 
  val = 0.0;
#endif

#if PROBLEM == 6 
  val = expProblemRhs( x ); 
#endif

#if PROBLEM == 7
  val = coscosHalf( x );
#endif

  val *= exactFactor();
  return val;
}

double exactPol(const double x[dim])
{
  double ret = 1.0;
  for(int i=0; i<dim; i++)
    ret *= ( x[i] - SQR(x[i]) );
  return ret;
}

double inSpringingCorner(const double x[dim])
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

double exactExpProblem(const double x[dim])
{
  double ret = 4.0 *(1.-SQR(x[0]))*(1.-SQR(x[1]))*exp(0.75*(x[0]+x[1]));
  return ret;
}

double exactSolution(const double x[dim])
{
  double val = 0.0;
#if PROBLEM == 2
#endif

#if PROBLEM == 3
  val = cos(M_PI * x[0]);
#endif

#if PROBLEM == 4 
  val = exactPol( x );
#endif

#if PROBLEM == 5 
  val = inSpringingCorner( x );
#endif

#if PROBLEM == 6 
#endif

#if PROBLEM == 7
  val = cos(0.5*M_PI*x[0]) * cos(0.5*M_PI*x[1]);
#endif

  val += globalShift;
  return val;
}


void exactPol(const double x[dim], double grad[dim])
{
  for(int j=0; j<dim; ++j) grad[j] = 1.0;

  for(int j=0; j<dim; ++j)
  {
    for(int i=0; i<dim; ++i)
    {
      if(i == j) grad[j] *= (1.0-2.0*x[i]);
      else       grad[j] *= ( x[i] - SQR(x[i]) );
    }
  }
}

void exactGradient(const double x[dim], double grad[dim])
{
#if PROBLEM == 2
#endif

#if PROBLEM == 3
  grad[0] = - M_PI * sin( M_PI * x[0] ); // -2.0;
  grad[1] = 0.0;
#endif

#if PROBLEM == 4 
  exactPol( x , grad );
#endif
  
#if PROBLEM == 5 
  grad[0] = grad[1] = 0.0; 
#endif
  
#if PROBLEM == 7
  grad[0] = 0.5*M_PI*-sin(0.5*M_PI*x[0])*cos(0.5*M_PI*x[1]);
  grad[1] = 0.5*M_PI*cos(0.5*M_PI*x[0])*-sin(0.5*M_PI*x[1]);
#endif

}

void rhsNeumann(const double x[dim], double grad[dim])
{
  exactGradient(x,grad);
  for(int i=0; i<dim; ++i) grad[i] *= exactFactor();
}

bool boundaryDataFunction(const double x[dim], double & val)
{
  abort();
#if PROBLEM == 3
  if(x[0] <= 0.0) 
  {
    val = globalShift + 1.0; 
    return true; 
  }

  if(x[0] >= 1.0)
  {
    val = globalShift - 1.0; 
    return true;
  }

  val = globalShift;
  return false;
#endif

#if PROBLEM == 4 
  val = globalShift;
  return true; 
#endif

#if PROBLEM == 5 
  val = inSpringingCorner( x ); 
  val += globalShift;
  return true; 
#endif

#if PROBLEM == 6 
  val = exactSolution( x ); 
  return true; 
#endif

#if PROBLEM == 7
  val = exactSolution( x ); 
  if(x[0] <= 0.0) return false;
  return true; 
#endif

  return true;
}
#endif

#endif
