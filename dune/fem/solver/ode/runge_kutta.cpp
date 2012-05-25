#include <cassert>
#include "ode_solver.hpp"
#include "blas.hpp"

using namespace pardg;



// class ExplicitRungeKutta
ExplicitRungeKutta::ExplicitRungeKutta(Communicator &comm, 
				       int num_of_stages, 
				       int order, 
				       Function &f) :
  ODESolver(comm, num_of_stages), f(f), num_of_stages(num_of_stages),
  A(num_of_stages), b(num_of_stages), c(num_of_stages), 
  alpha(num_of_stages-1), beta(num_of_stages), gamma(num_of_stages),
  order(order)
{}


ExplicitRungeKutta::ExplicitRungeKutta(Communicator &comm, 
				       int num_of_stages, 
				       int order, 
				       Function &f,
				       const double *a, 
				       const double *b, 
				       const double *c) :
  ODESolver(comm, num_of_stages), f(f), num_of_stages(num_of_stages),
  A(num_of_stages, a), b(num_of_stages, b), c(num_of_stages, c), 
  alpha(num_of_stages-1), beta(num_of_stages), gamma(num_of_stages),
  order(order)
{
  const int s = num_of_stages;

  // auxiliary matrices
  Matrix _A(s-1), _Ainv(s-1), _AL(s-1);
  for(int i=0; i<s-1; i++){
    for(int j=0; j<s-1; j++){
      _A(i,j) = A(i+1,j);
      if (i!=j) _AL(i,j) = _A(i,j);
    }
  }
  _Ainv = _A;
  _Ainv.inverse();

  // set alpha
  alpha = _AL * _Ainv;
  for(int i=0; i<s-1; i++) alpha(i,i) = _A(i,i);

  // set beta
  for(int i=0; i<s-1; i++){
    double sum = 0.0;
    for(int j=0; j<s-1; j++) sum += b[j] * _Ainv(j,i);
    beta[i] = sum;			       
  }
  beta[s-1] = b[s-1];

  // set gamma
  for(int i=0; i<s-1; i++){
    double sum = 1.0;
    for(int j=0; j<i; j++) sum -= alpha(i,j);
    gamma[i] = sum;
  }
  double sum = 1.0;
  for(int i=0; i<s-1; i++){
    for(int j=0; j<s-1; j++) sum -= b[i] * _Ainv(i,j);
  }
  gamma[s-1] = sum;

  if (false){
    std::cout << "alpha:" << std::endl << alpha << std::endl
	      << "beta:" << std::endl << beta << std::endl
	      << "gamma:" << std::endl << gamma << std::endl
	      << "c:" << std::endl << this->c << std::endl;
  }
}


bool ExplicitRungeKutta::step(double t, double dt, double *u)
{
  dim = f.dim_of_value();
  new_size(dim);
  return solve(t, dt, u);
}


inline
bool nearly_equal(double x, double y)
{
  const double epsilon = 1.0e-12;
  return fabs(x-y) < epsilon;
}

bool ExplicitRungeKutta::solve(double t, double dt, double *u)
{
  const int s = num_of_stages;
  const double *uim = u;
  double *ui = U;

  // build intermediate sulutions
  for(int i=0; i<s-1; i++){
    f(t + c[i]*dt, uim, ui);
    daxpby(dim, gamma[i], u, 1, dt*alpha(i,i), ui, 1);
    for(int j=0; j<i; j++){
      const double _alpha = alpha(i,j);
      if (!nearly_equal(_alpha, 0.0)){	
	cblas_daxpy(dim, _alpha, U+j*dim, 1, ui, 1);
      }
    }
    if (limiter) (*limiter)(ui);
    uim = ui;
    ui += dim;
  } 

  // update solution
  f(t + c[s-1]*dt, uim, ui);
  daxpby(dim, beta[s-1]*dt, ui, 1, gamma[s-1], u, 1);
  for(int i=0; i<s-1; i++) cblas_daxpy(dim, beta[i], U+i*dim, 1, u, 1);
  if (limiter) (*limiter)(u);

  return true;
}




// class ExplicitEuler
ExplicitEuler::ExplicitEuler(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 1, 1, f) {}


bool ExplicitEuler::solve(double t, double dt, double *u)
{
  f(t, u, U);
  cblas_daxpy(dim, dt, U, 1, u, 1);
  if (limiter != NULL) (*limiter)(u);

  return true;
}




// class ExplicitModifiedEuler
ExplicitModifiedEuler::ExplicitModifiedEuler(Communicator &comm,
					     Function &f) : 
  ExplicitRungeKutta(comm, 2, 2, f) {}


bool ExplicitModifiedEuler::solve(double t, double dt, double *u)
{
 f(t, u, U);
 daxpby(dim, 1.0, u, 1, 0.5*dt, U, 1);
 if (limiter) (*limiter)(U);

 double *U1 = U + dim;
 f(t + 0.5*dt, U, U1);
 cblas_daxpy(dim, dt, U1, 1, u, 1);
 if (limiter) (*limiter)(u);
 
 return true;
}




//class ExplicitTVD2, aka scheme by Heun
ExplicitTVD2::ExplicitTVD2(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 2, 2, f) {}


bool ExplicitTVD2::solve(double t, double dt, double *u)
{
 f(t, u, U);
 daxpby(dim, 1.0, u, 1, dt, U, 1);
 if (limiter) (*limiter)(U);

 double *U1 = U + dim;
 f(t + dt, U, U1);
 for(int k=0; k<dim; k++) u[k] = 0.5 * (u[k] + U[k] + dt*U1[k]);
 if (limiter) (*limiter)(u);
 
 return true;
}




//class ExplicitTVD3
ExplicitTVD3::ExplicitTVD3(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 3, 3, f) {}

bool ExplicitTVD3::solve(double t, double dt, double *u)
{
  f(t, u, U);
  daxpby(dim, 1.0, u, 1, dt, U, 1);
  if (limiter) (*limiter)(U);

  double *U1 = U + dim;
  f(t+dt, U, U1);
  for(int k=0; k<dim; k++) U1[k] = 0.75*u[k] + 0.25*(U[k] + dt*U1[k]);
  if (limiter) (*limiter)(U1);

  double *U2 = U1 + dim;
  f(t + 0.5*dt, U1, U2);
  for(int k=0; k<dim; k++){
    u[k] = (1.0/3.0) * ( u[k] + 2.0*(U1[k] + dt*U2[k]) );
  }
  if (limiter) (*limiter)(u);

  return true;
}

//class ExplicitRK3
static const double RK3_A[] =
  {0.0, 0.0, 0.0, 
   0.5, 0.0, 0.0, 
   -1.0, 2.0, 0.0, 
  };
static const double RK3_b[] = {1.0/6.0, 4.0/6.0, 1.0/6.0};
static const double RK3_c[] = {0.0, 0.5, 1.0};

ExplicitRK3::ExplicitRK3(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 3, 3, f, RK3_A, RK3_b, RK3_c) {}




//class ExplicitRK4b
static const double RK4_A[] =
  {0.0, 0.0, 0.0, 0.0,
   0.5, 0.0, 0.0, 0.0,
   0.0, 0.5, 0.0, 0.0,
   0.0, 0.0, 1.0, 0.0
  };
static const double RK4_b[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
static const double RK4_c[] = {0.0, 0.5, 0.5, 1.0};

ExplicitRK4b::ExplicitRK4b(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 4, 4, f, RK4_A, RK4_b, RK4_c) {}




 
//class ExplicitRK4
ExplicitRK4::ExplicitRK4(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 4, 4, f, RK4_A, RK4_b, RK4_c) {}


bool ExplicitRK4::solve(double t, double dt, double *u)
{
  const double dt_2 = 0.5 * dt;

  // u0 = u + dt/2 * f(u)
  f(t, u, U);
  daxpby(dim, 1.0, u, 1, dt_2, U, 1);
  if (limiter) (*limiter)(U);

  // u1 = u + dt/2 * f(u0)
  double *U1 = U + dim;
  f(t + 0.5*dt, U, U1);
  daxpby(dim, 1.0, u, 1, dt_2, U1, 1);
  if (limiter) (*limiter)(U1);

  // u2 = u + dt * f(u1)
  double *U2 = U1 + dim;
  f(t + 0.5*dt, U1, U2);
  daxpby(dim, 1.0, u, 1, dt, U2, 1);
  if (limiter) (*limiter)(U2);

  // u^{n+1} = 1/3 * ( -u^n + u0 + 2*u1 + u2 + 0.5*dt*f(u2) )
  double *U3 = U2 + dim;
  f(t + dt, U2, U3);
  for(int k=0; k<dim; k++){
    u[k] = (1.0/3.0) * ( U[k] + 2.0*U1[k] + U2[k] + dt_2*U3[k] - u[k] );
  }
  if (limiter) (*limiter)(u);

  return true;
}




//class ExplicitButcher6
static const double Butcher6_A[] =
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   2.0/9.0, 4.0/9.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   7.0/36.0, 2.0/9.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,
   -35.0/144.0, -55.0/36.0, 35.0/48.0, 15.0/8.0, 0.0, 0.0, 0.0,
   -1.0/360.0, -11.0/36.0, -1.0/8.0, 1.0/2.0, 1.0/10.0, 0.0, 0.0,
   -41.0/260.0, 22.0/13.0, 43.0/156.0, -118.0/39.0, 32.0/195.0, 80.0/39.0, 0.0
  };
static const double Butcher6_b[] = 
  {13.0/200.0, 0.0, 11.0/40.0, 11.0/40.0, 4.0/25.0, 4.0/25.0, 13.0/200.0};
static const double Butcher6_c[] = 
  {0.0, 0.5, 2.0/3.0, 1.0/3.0, 5.0/6.0, 1.0/6.0, 1.0};

ExplicitButcher6::ExplicitButcher6(Communicator &comm, Function &f) : 
  ExplicitRungeKutta(comm, 7, 6, f, Butcher6_A, Butcher6_b, Butcher6_c) {}


