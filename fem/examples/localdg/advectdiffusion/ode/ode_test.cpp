#include <iostream>
#include <cmath>
#include "ode_solver.hpp"
#include "function.hpp"
#include "linear_solver.hpp"

namespace DuneODE {
using namespace std;


double eoc(double h_fine, double h_coarse,
           double error_fine, double error_coarse)
{
  return log(error_fine/error_coarse) / log(h_fine/h_coarse);
}


// u'(t) = u(t)^2, general solution u(t) = -1/(t - c)
class SimplestFunction : public Function
{
public:
  SimplestFunction(int dim = 1) : dim(dim) {}

  void operator()(const double *u, double *f, int i = 0)
  {    
    if (i == 0){ // function
      for(int k=0; k<dim; k++) f[k] = u[k]*u[k];
      operation_count++;
    }
    else if (i == 1){ // first derivative
      for(int k=0; k<dim*dim; k++){
	if (k/dim == k%dim) f[k] = 1.0;
	else f[k] = 0.0;
      }
    }
    else assert(0);
  }

  int dim_of_argument(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int dim_of_value(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int operation_count;

private:
  const int dim;
};



// general solution of ODE: u(t) = c exp(t^2)
class SimpleFunction : public Function
{
public:
  SimpleFunction(int dim = 1) : dim(dim) {}

  void operator()(const double *u, double *f, int i = 0)
  {    
    if (i == 0){ // function
      for(int k=0; k<dim; k++) f[k] = 2.0 * time() * u[k];
      operation_count++;
    }
    else if (i == 1){ // first derivative
      for(int k=0; k<dim*dim; k++){
	if (k/dim == k%dim) f[k] = 2.0 * time();
	else f[k] = 0.0;
      }
    }
    else assert(0);
  }

  int dim_of_argument(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int dim_of_value(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int operation_count;

private:
  const int dim;
};



// SimpleFunction / frac
class FracSimpleFunction : public Function
{
public:
  FracSimpleFunction(double frac, int dim = 1) : frac(frac), dim(dim) {}

  void operator()(const double *u, double *f, int i = 0)
  {    
    if (i == 0){ // function
      for(int k=0; k<dim; k++) f[k] = frac * 2.0 * time() * u[k];
      operation_count++;
    }
    else if (i == 1){ // first derivative
      for(int k=0; k<dim*dim; k++){
	if (k/dim == k%dim) f[k] = frac * 2.0 * time();
	else f[k] = 0.0;
      }
    }
    else assert(0);
  }

  int dim_of_argument(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int dim_of_value(int i = 0) const 
  { 
    if (i==0) return dim;
    else if (i==1) return dim*dim;
    else assert(0);
  }

  int operation_count;

private:
  const double frac;
  const int dim;
};






int main(int argc, char *argv[])
{
  Communicator comm(argc, argv);

  const int dim = 1;
  SimpleFunction f(dim);
  const double frac = 0.3;
  FracSimpleFunction fim(frac, dim);
  FracSimpleFunction fex(1.0-frac, dim);
  double u[dim];

  ExplicitEuler ode_solver(comm, f);
  //ExplicitModifiedEuler ode_solver(comm, f);
  //ExplicitTVD2 ode_solver(comm, f);
  //ExplicitTVD3 ode_solver(comm, f);
  //ExplicitRK4 ode_solver(comm, f);
  //ExplicitRK4b ode_solver(comm, f);
  //ExplicitButcher6 ode_solver(comm, f);
  //ExplicitRK3 ode_solver(comm, f);
  //ExplicitSSP ode_solver(comm, f, 3);
  //DIRK3 ode_solver(comm, f);
  //ImplicitEuler ode_solver(comm, f);
  //SIRK23 ode_solver(comm, fim, fex);
  //SemiImplicitEuler ode_solver(comm, fim, fex);
  //IMEX_SSP222 ode_solver(comm, fim, fex);
  //SIRK33 ode_solver(comm, fim, fex);


  //ExplicitBulirschStoer ode_solver(comm, f, 5);
  //ImplicitBulirschStoer ode_solver(comm, f, 2);

  //LUSolver linear_solver;
  GMRES<20> linear_solver(comm);
  //linear_solver.DynamicalObject::set_output(cout);
  //linear_solver.IterativeSolver::set_output(cout);
  linear_solver.set_tolerance(1.0e-13);
  linear_solver.set_max_number_of_iterations(1000);

  //ode_solver.set_linear_solver(linear_solver);
  //ode_solver.set_tolerance(1.0e-8);
  ode_solver.DynamicalObject::set_output(cout);
  //ode_solver.IterativeSolver::set_output(cout);

  double error_old = 1.0;
  double dt = 0.01*1.0;
  const double sigma = 0.8;

  while (dt > 1.0e-3){
    dset(dim, 1.0, u, 1);
    //dset(dim, 0.5, u, 1);

    double T = 1.0;
    double t=0;
    f.operation_count = 0;

    while (t<T){
      const double t_new = t + dt;
      double dt_adj = dt;
      if (t_new > T) dt_adj = T-t;
      
      const bool convergence = ode_solver.step(t, dt_adj, u);
      assert(convergence);
      if (!convergence) break;
      
      if (t_new > T) t = T;
      else t = t_new;
    }

    double error = fabs(u[0] - exp(1.0));
    //double error = fabs(u[0] - 1.0);

    double order = eoc(dt, dt/sigma, error, error_old);    

    cout.width(14);
    cout << dt;
    cout.width(14);
    cout << error;
    cout.width(14);
    cout << order;
    cout.width(14);
    cout << f.operation_count << endl;

    error_old = error;
    dt *= sigma;
  }

  return 0;
}

}
