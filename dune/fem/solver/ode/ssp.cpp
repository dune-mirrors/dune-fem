#include "ode_solver.hpp"
#include "blas.hpp"


using namespace pardg;



static const double ssp_alpha[] =
  {1.0,
   0.5, 0.5,
   1.0/3.0, 1.0/2.0, 1.0/6.0,
   3.0/8.0, 1.0/3.0, 1.0/4.0, 1.0/24.0,
   11.0/30.0, 3.0/8.0, 1.0/6.0, 1.0/12.0, 1.0/120.0,
   53.0/144.0, 11.0/30.0, 3.0/16.0, 1.0/18.0, 1.0/48.0, 1.0/720.0
  };


ExplicitSSP::ExplicitSSP(Communicator &comm, Function &f, int m) :
  ODESolver(comm, 2), f(f), num_of_stages(m), 
  alpha(m, ssp_alpha+m*(m-1)/2)
{
  assert(m>0 && m<=6);
}



bool ExplicitSSP::step(double t, double dt, double *u)
{
  dim = f.dim_of_value();
  new_size(dim);
 
  double *ui = U;
  double *tmp = U + dim;
  cblas_dcopy(dim, u, 1, ui, 1);
  cblas_dscal(dim, alpha[0], u, 1);

  for(int i=0; i<num_of_stages-1; i++){
    f(t, ui, tmp);
    cblas_daxpy(dim, dt, tmp, 1, ui, 1); // u_i -> u_{i+1}
    cblas_daxpy(dim, alpha[i+1], ui, 1, u, 1);
  }
  f(t, ui, tmp);
  cblas_daxpy(dim, alpha[num_of_stages-1]*dt, tmp, 1, u, 1);

  return true;
}

