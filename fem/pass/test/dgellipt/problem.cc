#ifndef ELLIPTPROBLEM_CC
#define ELLIPTPROBLEM_CC

#include <math.h>

#define PROBLEM 1

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


double sinsin(const double x[dim]) 
{
  double sin_x = sin(2.0*M_PI*x[0]);
  double sin_y = sin(2.0*M_PI*x[1]);

  double val = -8.0 * M_PI * M_PI * sin_x * sin_y ;
  return -val;
}

double coscos(const double x[dim]) 
{
  double cos_x = cos(2.0*M_PI*x[0]);
  double cos_y = cos(2.0*M_PI*x[1]);
       
  double val = -8.0 * M_PI*M_PI* cos_x * cos_y ;
  return -val;
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
#if PROBLEM == 1 
  val = sinsin( x );
#endif

#if PROBLEM == 2
  val = coscos( x );
#endif

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
#if PROBLEM == 1 
  val = sin(2.0*M_PI*x[0]) * sin(2.0*M_PI*x[1]);
#endif

#if PROBLEM == 2
  val = cos(2.0*M_PI*x[0]) * cos(2.0*M_PI*x[1]);
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
  val = exactExpProblem( x );
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
#if PROBLEM == 1 
  grad[0] = 2.0*M_PI*cos(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1]);
  grad[1] = 2.0*M_PI*sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
#endif

#if PROBLEM == 2
  grad[0] = 2.0*M_PI*-sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
  grad[1] = 2.0*M_PI*cos(2.0*M_PI*x[0])*-sin(2.0*M_PI*x[1]);
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
#if PROBLEM == 1 
  val = exactSolution( x ); 
  if(x[0] <= 0.0) return false;
  return true; 
#endif

#if PROBLEM == 2
  val = exactSolution( x ); 
  if(x[0] <= 0.0) return false;
  return true; 
#endif

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
