/*******************************************************************
 Problem description in domain D_gdl + D_cat:

    M4) Mass and momentum balance for liquid and gas flow in a porous media
           -input   variables:
           -free    variables: water saturation (s_w), global pressure (p), global velocity (u)
           -derived variables: 

  ****************************************************************
  
*/
#ifndef FUELCELL_FCM4_CC
#define FUELCELL_FCM4_CC

#include <math.h>

#define PROBLEM 1

// shift solutions to avoid dirichlet 0 bnd 
static const double globalShift = 2.0;

double exactFactor() 
{
  return 1.0/(5.24*1e5);
  //return 2.0;
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
  for(int i=0; i<dim; ++i) grad[i] *= exactFactor();
}

void rhsNeumann(const double x[dim], double grad[dim])
{
  exactGradient(x,grad);
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
  //if(x[0] <= 0.0) return false;
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


}

namespace Twophase {


  

//! the residual saturations 
static const double residual_saturation_w = 0.2; 
static const double residual_saturation_g = 0.15; 

static const double residual_saturation_w_g_1 = 1.0/0.65;

// relative permeabilty kr_i , i=w,g  
/* kr_w/g(s)           : relative permeability of water/gas */

/* mue_w/g                : viscosity of water (1.002e-3 Pa*s at 20 deg. Celc.) 
                                         gas   (1.720e-5 Pa*s at  0 deg. Celc.) */
static const double mue_w   = 0.001; 
static const double mue_w_1 = 1.0/mue_w; 
static const double mue_g   = 0.01; // here the non wetting phase 
static const double mue_g_1 = 1.0/mue_g;

double residual_saturation(double s)
{
  assert( s >= 0.0 );
  assert( s <= 1.0 );
  return (s - residual_saturation_w) * residual_saturation_w_g_1;
}

//! realive permeability of water phase (van Genuchten)
//tex: \newcommand{\krwgdl}{$kr_w (s) = (s_e(s))^4} 
//tex: The relative permeabilty of water is defined as \krwgdl.
double kr_w (double u)
{
  double s = residual_saturation(u);
  return (SQR(s)*SQR(s));
}

//! realtive permeabilty of gas phase (van Genuchten)
//tex: The relative permeabilty of gas is defined $kr_g(s) = kr_w(1-s)$.
double kr_g (double u)
{
  double s = residual_saturation(u);
  return SQR(1.0 - s) * (1.0 - SQR(s));
}

/*! lambda_w (s)  = kr_w(s) / mue_w    : relative mobility of water */
//tex: The relative mobiltyy of water is defined $\lambda_w(s) = kr_w(s)/\mue_w$. 
double lambda_w (double s)
{
  return kr_w(s) * mue_w_1;
}

double DF_lambda_w (double u)
{
  double s = residual_saturation(u);
  return ((mue_w_1) * 4.0 * SQR(s)*s) * residual_saturation_w_g_1;
}

/*! lambda_g (s)  = kr_g(s) / mue_g    : relative mobility of gas */
//tex: The relative mobiltyy of water is defined $\lambda_g(s) = kr_g(s)/\mue_g$. 
double lambda_g (double s)
{
  return kr_g(s) * mue_g_1;
}

double DF_lambda_g (double u)
{
  double s = residual_saturation(u);
  return (-mue_g_1) * ((2.0*(1.0-s)*(1.0-SQR(s))) + (SQR(1.0-s)*(2.0*s))) * residual_saturation_w_g_1;
}

/*! lambda (s)  = lambda_w(s) + lambda_g(s) : total mobility */
//tex: The total mobilty is defined as the sum of the relative mobilty, $\lambda(s) = \lambda_w(s) + \lambda_g(s)$. 
double lambda (double s)
{
  return lambda_w(s) + lambda_g(s);
}

/*! f_w (s)  = lambda_w(s) / lambda_w(s) : fractional flow rate of water */
//tex: The fractional flow rate of water is $f_w(s) = \frac{\lambda_w(s)}{\lambda(s)}$. 
double f_w (double s)
{
  return lambda_w(s) / lambda(s);
}

double Df_w (double s)
{
  // the residual saturation is evaluated in each function 
  double dfw = DF_lambda_w(s);

  double dfg = DF_lambda_g(s);
  double w   = lambda_w (s);
  double g   = lambda_g (s);

  double sumwn_2 = SQR(w + g);
  return (dfw*g - w*dfg)/sumwn_2;
}

//! porosity 
double porosity (const double x[dim], const double time )
{
  return 0.2;
}

static double velocityValue = 9e-7 / 0.2; // velo/porosity 

//! constant velocity 
void velocity (double velo[dim])
{
  velo[0] = velocityValue; // == 9e-7 / 0.2 
  for(int i=1; i<dim; i++) velo[i] = 0.0;
  return;
}      

//***********************************************************************

static const double permeabilty = 1e-11; 
void K (const double x[dim], const double time, double res [dim][dim] )
{
  for(int i=0 ;i<dim ; i++) 
    for(int j=0; j<dim ; j++) 
    {
      res[i][j] = 0.0; 
      if( i == j ) res[i][j] = permeabilty; 
    }
}

void k (const double x[dim], const double time, double res [1] )
{
  res[0] = permeabilty;
}

/* Parameter in the Brooks-Corey model for the capillary pressure */
/* see paper Peter Bastian  */
static const double BrooksCorey_pd = 1e4; // entry pressure 
static const double BrooksCorey_lambda = -0.5; // parameter lambda 

//tex: The capilary pressure--saturation relation ship is defined as 
//tex: $p_c(s) = -p_d (1-s)^{\frac{-1}{l_a}}$. 
double capillary_pressure (const double s)
{
  return BrooksCorey_pd  * pow( s , BrooksCorey_lambda );
}

double initial_saturation(const double x[dim])
{
  //if(x[0] < 150.0) return 1.0 - residual_saturation_g;
  return residual_saturation_w;
}

double dirichlet_saturation_left(const double x[dim])
{
  return 1.0 - residual_saturation_g;
}

double dirichlet_pressure_upper(const double x[dim])
{
  return x[0] * velocityValue;
  //return 3.2e5;
}

double dirichlet_pressure_lower(const double x[dim])
{
  return x[0] * velocityValue;
  //return 1.0e5;
}

double dirichlet_saturation_right(const double x[dim])
{
  return residual_saturation_w;
}

} // end namespace FuelCell

#endif
