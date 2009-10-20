#include <cassert>
#include <cfloat>
#include "ode_solver.hpp"
#include "blas.hpp"


using namespace pardg;



// class ExplicitBulirschStoer

// F = {2, 4, 8, 16, 32, ...}, F(0)=2
int ExplicitBulirschStoer::DoubleRombergSequence(int i)
{
  return 2 << (i+1);
}


// F = {2, 4, 6, 8, 10, ...}, F(0)=2
int ExplicitBulirschStoer::DoubleHarmonicSequence(int i)
{
  return 2*(i+1);
}


ExplicitBulirschStoer::ExplicitBulirschStoer(Communicator &comm, 
					     Function &f, 
					     int num_of_stages, 
					     int (*seq)(int) ) :
  ODESolver(comm, num_of_stages+2), num_of_stages(num_of_stages), 
  f(f), sequence(seq)
{}


// Strehmel, Karl;  Weiner, Rüdiger
// Numerik gewöhnlicher Differentialgleichungen. 
// (Numerical methods for ordinary differential equations). (German)
// Teubner Studienbücher: Mathematik. Stuttgart: Teubner. 462 p. (1995). 
// S. 77 ff
//
// Stoer, Josef;  Bulirsch, Roland
// Numerische Mathematik II. Eine Einführung - unter Berücksichtigung 
// von Vorlesungen von F. L. Bauer. 3., verb. Aufl. (German)
// Springer-Lehrbuch. Berlin etc.: Springer-Verlag. xiii, 341 p. (1990).
// S. 166
bool ExplicitBulirschStoer::step(double t, double dt, double *u)
{
  dim = f.dim_of_value();
  new_size(dim);
  double *tmp = U + (num_of_stages+1)*dim;

  // store f(t,u) in U[0], this is necessary in every stage,
  // can be (and will be) overwritten in the last stage
  f(t, u, U); 

  for(int i=0; i<num_of_stages; i++){
    const int n = sequence(i); // n_i

    // compute u_{dt/n_i} (t+dt) by application of the modified midpoint 
    // scheme, store it in U[num_of_stages-1-i]

    // init
    const double dt_n = dt/n;
    double *z_km = U + (num_of_stages-1-i)*dim;
    double *z_k = U + num_of_stages*dim;
    dwaxpby(dim, 1.0, u, 1, dt_n, U, 1, z_k, 1); // z_k = u + dt/n * f(t,u)
    if (limiter) (*limiter)(z_k);
    cblas_dcopy(dim, u, 1, z_km, 1);
  
    // loop
    for(int k=1; k<n; k++){
      f(t + k*dt_n, z_k, tmp);
      double *z_kp = z_km;
      cblas_daxpy(dim, 2.0*dt_n, tmp, 1, z_kp, 1);
      if (limiter) (*limiter)(z_kp);
      z_km = z_k;
      z_k = z_kp;
    }

    // finalize
    f(t + dt, z_k, tmp);
    double *_U = U + (num_of_stages-1-i)*dim;
    for(int l=0; l<dim; l++) _U[l] = 0.5*( z_km[l] + z_k[l] + dt_n*tmp[l] );
    if (limiter) (*limiter)(_U);


    // Extrapolation process
    // T_{i,0} = u_{dt/n_i} is stored in U[num_of_stages-1-i],
    // compute T_{i,1},..,T_{i,i} stored in 
    // U[num_of_stages-1-i+1],..,U[num_of_stages-1]

    for(int j=1; j<=i; j++){
      const double h_i_j = 1.0 / sequence(i-j); // 1.0/n_{i-j}      
      const double alpha = 1.0 / (n*n*h_i_j*h_i_j - 1.0);
      daxpby(dim, 1.0+alpha, U+(num_of_stages-1-i+j-1)*dim, 1, 
	     -alpha, U+(num_of_stages-1-i+j)*dim, 1);
    }    
  }

  // update approximate solution, stored in U[num_of_stages-1] after 
  // num_of_stages stages
  cblas_dcopy(dim, U+(num_of_stages-1)*dim, 1, u, 1);

  return true;
}
  



// class ImplicitBulirschStoer

int ImplicitBulirschStoer::DoubleRombergSequence(int i)
{
  return ExplicitBulirschStoer::DoubleRombergSequence(i);
}


int ImplicitBulirschStoer::DoubleHarmonicSequence(int i)
{
  return ExplicitBulirschStoer::DoubleHarmonicSequence(i);
}


void ImplicitBulirschStoer::set_linear_solver(DirectLinearSolver &ls)
{
  dls = &ls;
}


void ImplicitBulirschStoer::set_linear_solver(IterativeLinearSolver &ls)
{
  ils = &ls;
}



void ImplicitBulirschStoer::resize(int new_size, int component)
{
  // new_size >= dim
  delete[] U;
  if (dls){ // for direct linear solver use
    U = new double[ (num_of_stages + 3 + new_size)*new_size ]; 
  }
  else if (ils){
    U = new double[ (num_of_stages + 4)*new_size];
  }

  tmp = U + (num_of_stages+1)*new_size;
  F = tmp + new_size;  
}




ImplicitBulirschStoer::ImplicitBulirschStoer(Communicator &comm, 
					     Function &f, 
					     int num_of_stages, 
					     int (*seq)(int) ) :
  ODESolver(comm, 0), num_of_stages(num_of_stages), 
  f(f), sequence(seq), dls(NULL), ils(NULL), op(*this)
{
  // set this to some useful values
  tolerance = 1.0e-6;
  max_num_of_iterations = 20;
}



bool ImplicitBulirschStoer::step(double t, double dt, double *u)
{
  dim = f.dim_of_value();
  new_size(dim);

  if (dls){ // use direct linear solver, serial only
    assert( comm.size() == 1 );
    return step_direct(t, dt, u);
  }
  else if (ils){ // use iterative linear solver
    return step_iterative(t, dt, u);
  }
  else assert(0); // linear solver needed

  return false;
}




bool ImplicitBulirschStoer::step_iterative(double t, double dt, double *u)
{
  double *p = F + dim;
  this->t = t;

  for(int i=0; i<num_of_stages; i++){
    const int n = sequence(i); // n_i

    // compute u_{dt/n_i} (t+dt) by application of the implicit midpoint
    // scheme, store it in U[num_of_stages-1-i]

    // init
    dt_n = dt/n;
    const double dt_2n = 0.5*dt_n;
    z_km = U + (num_of_stages-1-i)*dim;
    z_k = U + num_of_stages*dim;
    cblas_dcopy(dim, u, 1, z_km, 1); // z_km = u^n

    // loop
    for(k=0; k<n; k++){
    
      if (true) cblas_dcopy(dim, z_km, 1, z_k, 1); // z_k = z_km
      else { // z_k = z_km + dt/n f(t+(k+0.5)*dt/n, z_k) euler pedictor step
	f(t + (k+0.5)*dt_n, z_k, tmp);
	dwaxpby(dim, 1.0, z_km, 1, dt_2n, tmp, 1, z_k, 1);
      }


      // Newton iteration for solving 
      // z_k = z_km + dt/n * f(t + (k+0.5)*dt/n, 0.5*(z_k + z_km))
      int num_of_iterations = 0;
      double dist_old = DBL_MAX;
      while (true) {
	for(int l=0; l<dim; l++) tmp[l] = 0.5 * (z_k[l] + z_km[l]); 
	f(t + (k+0.5)*dt_n, tmp, F);
	for(int l=0; l<dim; l++) F[l] = z_k[l] - z_km[l] - dt_n*F[l];

	dset(dim, 0.0, p, 1);
	bool lin_solver_conv = ils->solve(op, p, F);
	if (!lin_solver_conv) return false;

	double global_dot, local_dot = 0.0;
	for(int l=0; l<dim; l++){
	  z_k[l] -= p[l];
	  local_dot += p[l]*p[l];
	}
	comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
	const double dist = sqrt(global_dot);

	if (IterativeSolver::os){
	  *IterativeSolver::os << "Newton: iteration: "
			       << num_of_iterations << "    "
			       << "|p|: " << dist << "   "
			       << std::endl;
	}
	
	if (dist < tolerance){ // successful solving 
	  double *z_kp = z_km;
	  z_km = z_k;
	  z_k = z_kp;
	  break;
	}
	else if (num_of_iterations >= max_num_of_iterations 
		 || dist >= dist_old) return false; // not successful
	dist_old = dist;
	num_of_iterations++;
      }
      // end of Newton iteration, approx solution is stored in z_km

      // store approx solution at (k+1)*dt/n in the right position
      // if it isnt (because of swapping of variables z_k and z_km)
      double *_U = U + (num_of_stages-1-i)*dim;
      if (z_km != _U) cblas_dcopy(dim, z_km, 1, _U, 1);
    }


    // Extrapolation process
    // T_{i,0} = u_{dt/n_i} is stored in U[num_of_stages-1-i],
    // compute T_{i,1},..,T_{i,i} stored in
    // U[num_of_stages-1-i+1],..,U[num_of_stages-1]

    for(int j=1; j<=i; j++){
      const double h_i_j = 1.0 / sequence(i-j); // 1.0/n_{i-j}
      const double alpha = 1.0 / (n*n*h_i_j*h_i_j - 1.0);
      daxpby(dim, 1.0+alpha, U+(num_of_stages-1-i+j-1)*dim, 1,
             -alpha, U+(num_of_stages-1-i+j)*dim, 1);
    }
  }


  // update approximate solution, stored in U[num_of_stages-1] after
  // num_of_stages stages
  cblas_dcopy(dim, U+(num_of_stages-1)*dim, 1, u, 1);

  return true;
}
  


bool ImplicitBulirschStoer::step_direct(double t, double dt, double *u)
{
  double *DF = F + dim;

  for(int i=0; i<num_of_stages; i++){
    const int n = sequence(i); // n_i

    // compute u_{dt/n_i} (t+dt) by application of the implicit midpoint
    // scheme, store it in U[num_of_stages-1-i]

    // init
    dt_n = dt/n;
    const double dt_2n = 0.5*dt_n;
    z_km = U + (num_of_stages-1-i)*dim;
    z_k = U + num_of_stages*dim;
    cblas_dcopy(dim, u, 1, z_km, 1); // z_km = u^n

    // loop
    for(k=0; k<n; k++){
    
      if (true) cblas_dcopy(dim, z_km, 1, z_k, 1); // z_k = z_km
      else { // z_k = z_km + dt/n f(t+(k+0.5)*dt/n, z_k) euler pedictor step
	f(t + (k+0.5)*dt_n, z_k, tmp);
	dwaxpby(dim, 1.0, z_km, 1, dt_2n, tmp, 1, z_k, 1);
      }


      // Newton iteration for solving 
      // z_k = z_km + dt/n * f(t + (k+0.5)*dt/n, 0.5*(z_k + z_km))
      int num_of_iterations = 0;
      double dist_old = DBL_MAX;
      while (true) {
	for(int l=0; l<dim; l++) tmp[l] = 0.5 * (z_k[l] + z_km[l]); 
	f(t + (k+0.5)*dt_n, tmp, F);
	for(int l=0; l<dim; l++) F[l] = z_k[l] - z_km[l] - dt_n*F[l];

	f(t + (k+0.5)*dt_2n, tmp, DF, 1);
	for(int l=0; l<dim*dim; l++){
	  if (l/dim == l%dim) DF[l] = (1.0 - dt_2n*DF[l]);
	  else DF[l] *= -dt_2n;
	}

	dls->prepare(dim, DF);
	dls->solve(F); // DF dz = F, dz stored in F

	double dist = 0.0;
	for(int l=0; l<dim; l++){
	  z_k[l] -= F[l];
	  dist += F[l]*F[l];
	}
	dist = sqrt(dist);

	if (IterativeSolver::os){
	  *IterativeSolver::os << "Newton: iteration: "
			       << num_of_iterations << "    "
			       << "|du|: " << dist << "   "
			       << std::endl;
	}
	
	if (dist < tolerance){ // successful solving 
	  double *z_kp = z_km;
	  z_km = z_k;
	  z_k = z_kp;
	  break;
	}
	else if (num_of_iterations >= max_num_of_iterations 
		 || dist >= dist_old) return false; // not successful
	dist_old = dist;
	num_of_iterations++;
      }
      // end of Newton iteration, approx solution is stored in z_km

      // store approx solution at (k+1)*dt/n in the right position
      // if it isnt (because of swapping of variables z_k and z_km)
      double *_U = U + (num_of_stages-1-i)*dim;
      if (z_km != _U) cblas_dcopy(dim, z_km, 1, _U, 1);
    }


    // Extrapolation process
    // T_{i,0} = u_{dt/n_i} is stored in U[num_of_stages-1-i],
    // compute T_{i,1},..,T_{i,i} stored in
    // U[num_of_stages-1-i+1],..,U[num_of_stages-1]

    for(int j=1; j<=i; j++){
      const double h_i_j = 1.0 / sequence(i-j); // 1.0/n_{i-j}
      const double alpha = 1.0 / (n*n*h_i_j*h_i_j - 1.0);
      daxpby(dim, 1.0+alpha, U+(num_of_stages-1-i+j-1)*dim, 1,
             -alpha, U+(num_of_stages-1-i+j)*dim, 1);
    }
  }


  // update approximate solution, stored in U[num_of_stages-1] after
  // num_of_stages stages
  cblas_dcopy(dim, U+(num_of_stages-1)*dim, 1, u, 1);

  return true;
}
  



// ======== ImplicitBulirschStoer::LinearOperator implementation
ImplicitBulirschStoer::LinearOperator::
LinearOperator(ImplicitBulirschStoer &ibs) : ibs(ibs) {}


int ImplicitBulirschStoer::LinearOperator::dim_of_argument(int i) const
{
  return ibs.dim;
}


int ImplicitBulirschStoer::LinearOperator::dim_of_value(int i) const
{
  return ibs.dim;
}


void ImplicitBulirschStoer::LinearOperator::operator()(const double *p,
						       double *DFu_p, int i)
{
  const int dim = ibs.dim;
  const double *u = ibs.z_k;
  double local_dot[2], global_dot[2];
  local_dot[0] = cblas_ddot(dim, u, 1, u, 1);
  local_dot[1] = cblas_ddot(dim, p, 1, p, 1);
  ibs.comm.allreduce(2, local_dot, global_dot, MPI_SUM);
  const double norm_u = sqrt(global_dot[0]);
  const double norm_p_sq = global_dot[1];

  const double eps = (norm_p_sq > DBL_EPSILON)?
    sqrt( (1.0+norm_u)*DBL_EPSILON / norm_p_sq ) : sqrt(DBL_EPSILON);
  const double eps_inv = 1.0 / eps;

  // DFu_p = ( F(u+eps*p) - F(u) ) / eps
  // whith F(u) = u - z_km - dt/n * f( t + (k+0.5)*dt/n, 0.5*(u + z_km) )
  for(int i=0; i<dim; i++) ibs.tmp[i] = 0.5*(u[i] + eps*p[i] + ibs.z_km[i]);
  (ibs.f)(ibs.t + (ibs.k+0.5)*ibs.dt_n, ibs.tmp, DFu_p);
  for(int i=0; i<dim; i++){
    DFu_p[i] = eps_inv*(u[i]+eps*p[i]-ibs.z_km[i]-ibs.dt_n*DFu_p[i]-ibs.F[i]);
  }
}
