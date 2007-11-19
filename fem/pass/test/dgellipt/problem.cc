#ifndef ELLIPTPROBLEM_CC
#define ELLIPTPROBLEM_CC

#include <cmath>
#include <cassert>

template<int dim, class DomainField, class Field> 
class DataFunctionIF
{
protected:  
  DataFunctionIF() {}
public:  
  virtual ~DataFunctionIF() {}
  // diffusion tensor 
  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const = 0;
  // right hand side 
  virtual Field rhs  (const DomainField x[dim]) const = 0;
  // exact solution 
  virtual Field exact(const DomainField x[dim]) const = 0;
  // exact gradient 
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const = 0;

  // boundary data 
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const = 0;

  // neumann boundary 
  virtual void neumann(const DomainField x[dim], Field grad[dim]) const 
  {
    Field tmp[dim];
    gradExact(x,tmp);
    Field k[dim][dim];
    K(x,k);
    for(int i=0; i<dim; ++i) 
    {
      grad[i] = 0;
      for(int j=0; j<dim; ++j) 
        grad[i] += tmp[j] * k[i][j];
    }
  }
};


template <int dim, class DomainField, class Field> 
class SinSin : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
public:  
  virtual ~SinSin() {}
  SinSin(Field globalShift, Field factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    //assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }

  virtual Field rhs  (const DomainField x[dim]) const 
  {
    Field sin_x = sin(2.0*M_PI*x[0]);
    Field sin_y = sin(2.0*M_PI*x[1]);

    Field val = 8.0 * M_PI * M_PI * sin_x * sin_y ;
    val *= factor_;
    return val;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = sin(2.0*M_PI*x[0]) * sin(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*cos(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    if(x[0] <= 0.0) return false;
    return true; 
  }
};

template <int dim, class DomainField, class Field> 
class CosCos : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
public:  
  virtual ~CosCos() {}
  CosCos(Field globalShift, Field factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    //assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }


  virtual Field rhs  (const DomainField x[dim]) const 
  {
    Field cos_x = cos(2.0*M_PI*x[0]);
    Field cos_y = cos(2.0*M_PI*x[1]);
       
    Field val = 8.0 * M_PI*M_PI* cos_x * cos_y ;
    val *= factor_;
    return val;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = cos(2.0*M_PI*x[0]) * cos(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*-sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*cos(2.0*M_PI*x[0])*-sin(2.0*M_PI*x[1]);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    if(x[0] <= 0.0) return false;
    return true; 
  }
};

//! problem from Castillo paper 
template <int dim, class DomainField, class Field> 
class CastilloProblem : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
public:  
  virtual ~CastilloProblem() {}
  CastilloProblem(Field globalShift, Field factor)
    : globalShift_(globalShift)
    , factor_(factor) 
  {
    assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }


  virtual Field rhs  (const DomainField x[dim]) const 
  {
    Field ret = 0.0;
    Field tmp =-23+7*SQR(x[1])-24*x[0]+24*x[0]*SQR(x[1])+7*SQR(x[0])+9*SQR(x[0])*SQR(x[1])-24
                *x[1]+24*x[1]*SQR(x[0]);

    ret=-0.5*exp(0.75*(x[0]+x[1]));
    ret *= tmp;
    ret *= factor_;
    return ret;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = 4.0 *(1.-SQR(x[0]))*(1.-SQR(x[1]))*exp(0.75*(x[0]+x[1]));
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    // not implemented yet 
    grad[0] = grad[1] = 0.0; 
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};

//! problem Einstringende Ecke 
template <int dim, class DomainField, class Field> 
class InSpringingCorner : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
public:  
  virtual ~InSpringingCorner() {}
  InSpringingCorner(Field globalShift, Field factor)
    : globalShift_(0.0)
    , factor_(1.0) 
  {
    assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }

  virtual Field rhs  (const DomainField x[dim]) const 
  {
    return 0.0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field ret = 0.0;
    const Field Pi= M_PI; 
    Field tmp2=0.0;

    Field s = 0.0 ;
    for(int i=0; i<dim; ++i) s += x[i]*x[i];
    Field r = sqrt(s); 
   
    Field phi=0.0;
    Field fac=2.0;
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
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    // not implemented yet 
    grad[0] = grad[1] = 0.0; 
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


//! problem Riviere-Bastian 
template <int dim, class DomainField, class Field> 
class RiviereProblem : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
public:  
  virtual ~RiviereProblem() {}
  RiviereProblem(Field globalShift, Field factor)
    : globalShift_(0.0)
    , factor_(1.0) 
  {
    assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }

  virtual Field rhs  (const DomainField x[dim]) const 
  {
    Field val = exact(x);
    Field x_part = val * SQR(-2.0 * (x[0] - 0.5)) - 2.0 * val; 
    Field y_part = val * SQR(-2.0 * (x[1] - 0.5)) - 2.0 * val; 
    return -(x_part + y_part);
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field power = -( SQR( x[0] - 0.5 ) + SQR( x[1] - 0.5 ) );
    return pow( M_E , power );
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    Field val = exact( x );
    grad[0] = val * ( -2.0*x[0] + 1.0 );
    grad[1] = val * ( -2.0*x[1] + 1.0 );
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};

template <int dim, class DomainField, class Field> 
class CompactDG : public DataFunctionIF<dim,DomainField,Field>
{
  const Field globalShift_;
  const Field factor_;
  const Field alpha_;
  const Field beta_;
  const Field a_;
  const Field b_;
  const Field c_;
  const Field d_;
public:  
  virtual ~CompactDG() {}
  CompactDG(Field globalShift, Field factor)
    : globalShift_(globalShift)
    , factor_(factor) 
    , alpha_( 0.1 ) 
    , beta_ ( 0.3 )
    , a_ ( 5.1 ), b_ (-6.2 ), c_(4.3), d_(3.4)
  {
    //assert(dim == 2);
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const 
  {
    for(int i=0; i<dim; ++i) 
    {
      k[i][i] = factor_;
      for(int j=0; j<i; ++j)     k[i][j] = 0;
      for(int j=i+1; j<dim; ++j) k[i][j] = 0;
    }
  }


  virtual Field rhs  (const DomainField x[dim]) const 
  {
    Field val = 0.0;
    val *= factor_;
    return val;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = exp( alpha_ * sin( a_ * x[0] + b_ * x[1] ) + beta_ * cos( c_ * x[0] + d_ * x[1] ));
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 0.0; 
    grad[1] = 0.0; 
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


#endif
