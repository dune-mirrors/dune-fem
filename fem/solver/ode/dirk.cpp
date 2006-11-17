#include <cfloat>
#include "ode_solver.hpp"


using namespace pardg;



DIRK::DIRK(Communicator &comm, int num_of_stages, int order, Function &f,
	   const double *a, const double *b, const double *c) :
  ODESolver(comm, 0), f(f), num_of_stages(num_of_stages),
  order(order),
  A(num_of_stages, a), b(num_of_stages, b), c(num_of_stages, c),
  alpha(num_of_stages), beta(num_of_stages), gamma(num_of_stages),
  F(NULL), y(NULL), 
  ils(NULL), op(comm, f, dim, u_tmp, f_tmp), u_tmp(NULL), f_tmp(NULL)
{
  // set this to some useful values
  tolerance = 1.0e-6;
  max_num_of_iterations = 20;

  // build matrix alpha
  Matrix Ainv = A;
  Ainv.inverse();
  Matrix AL = A;
  for(int i=0; i<num_of_stages; i++) AL(i,i) = 0.0;
  alpha = AL * Ainv;
  for(int i=0; i<num_of_stages; i++) alpha(i,i) = A(i,i);
  
  // build vector beta
  for(int i=0; i<num_of_stages; i++){
    beta[i] = 0.0;
    for(int j=0; j<num_of_stages; j++) beta[i] += b[j] * Ainv(j,i);
  }
  
  // build vector gamma
  for(int i=0; i<num_of_stages; i++){
    gamma[i] = 1.0;
    for(int j=0; j<i; j++) gamma[i] -= alpha(i,j);
  }

  // delta
  delta = 1.0;
  for(int i=0; i<num_of_stages; i++) delta -= beta[i];
}


void DIRK::set_linear_solver(IterativeLinearSolver &ls)
{
  ils = &ls;
}


void DIRK::resize(int new_size, int component)
{
  // new_size >= dim
  delete[] U;
  U = new double[ (num_of_stages + 5) * new_size ];
  assert(U);

  dset((num_of_stages + 5) * new_size, 0.0, U, 1);

  Fpre = U + num_of_stages * new_size;
  F = Fpre + new_size;    
  y = F + new_size;
  u_tmp = y + new_size;    // for LinearOperator
  f_tmp = u_tmp + new_size; // 
}



bool DIRK::step(double t, double dt, double *u)
{
  dim = f.dim_of_value();
  new_size(dim);

  const bool convergence =  step_iterative(t, dt, u);

  // update solution
  if (convergence){
    cblas_dscal(dim, delta, u, 1);
    for(int i=0; i<num_of_stages; i++){
      cblas_daxpy(dim, beta[i], U+i*dim, 1, u, 1);
    }
    return true;
  }
  else return false;
}


bool DIRK::step_iterative(double t, double dt, double *u)
{
  for(int i=0; i<num_of_stages; i++){
    double *ui = U+i*dim;

    // prediction uf ui, todo: extrapolation or something...
    // Euler predictor
    f(t+c[i]*dt, u, f_tmp);
    dwaxpby(dim, c[i]*dt, f_tmp, 1, 1.0, u, 1, ui, 1);    
    // ui = u^n
    //cblas_dcopy(dim, u, 1, ui, 1);

    // setup of Fpre
    const double _gamma = gamma[i];
    for(int l=0; l<dim; l++) Fpre[l] = _gamma * u[l];
    for(int j=0; j<i; j++) cblas_daxpy(dim, alpha(i,j), U+j*dim, 1, Fpre, 1);

    // Newton iteration
    int iterations = 0;
    while (iterations < max_num_of_iterations){
      // setup f_tmp and F
      f(t + c[i]*dt, ui, f_tmp);
      const double lambda = alpha(i,i) * dt;
      for(int l=0; l<dim; l++) F[l] = ui[l] - lambda*f_tmp[l] - Fpre[l];

      // solve linear system
      dset(dim, 0.0, y, 1);
      op.setup(t+c[i]*dt, ui, lambda);
      const bool lin_solver_conv = ils->solve(op, y, F);
      if (!lin_solver_conv) return false;

      // update ui & apply limiter
      cblas_daxpy(dim, -1.0, y, 1, ui, 1);
      if (limiter) (*limiter)(ui);

      // compute norm_y
      double global_dot, local_dot;
      local_dot = cblas_ddot(dim, y, 1, y, 1);
      comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);

      if (IterativeSolver::os){
	*IterativeSolver::os << "Newton: iteration: "
			     << iterations << "    "
			     << "|p|: " << sqrt(global_dot) << "   "
			     << std::endl;
      }

      iterations++;    

      if(sqrt(global_dot) < tolerance) break;      
    }

    if (iterations >= max_num_of_iterations) return false;    
  }

  return true;
}


// bool DIRK::step_iterative(double t, double dt, double *u)
// {
//   // set f_tu
//   f(t, u, f_tu);

//   const int s = num_of_stages;
//   int iterations = 0;

//   // Newton iteration
//   while (iterations < max_num_of_iterations){
//     // setup of F
//     //std::cout << "U[0]: " << U[0] << std::endl;
//     for(int i=0; i<s; i++){
//       dwaxpby(dim, 1.0, U+i*dim, 1, -1.0, u, 1, F+i*dim, 1);
//     }
//     for(int j=0; j<s; j++){
//       f(t + c[j]*dt, U+j*dim, y);
//       for(int i=0; i<s; i++) cblas_daxpy(dim, -dt*A(i,j), y, 1, F+i*dim, 1);
//     }

//     // solve linear system & update U
//     double sq_norm_y = 0;
//     for(int i=0; i<s; i++){
//       // solve the linear system
//       const double *Fi = F + i*dim;
//       dset(dim, 0.0, y, 1);
//       op.setup(t, u, dt*A(i,i));
//       const bool lin_solver_conv = ils->solve(op, y, Fi);
//       //std::cout << lin_solver_conv << std::endl;
//       //for(int l=0; l<dim; l++) std::cout << y[l] << std::endl;
//       if (!lin_solver_conv) return false;

//       // update Fk
//       for(int k=i+1; k<s; k++){
// 	double *Fk = F + k*dim;
// 	const double lambda = A(i,k)/A(i,i);
// 	for(int l=0; l<dim; l++) Fk[l] -= lambda * (y[l] - Fi[l]);
//       }

//       // update ui
//       cblas_daxpy(dim, -1.0, y, 1, U+i*dim, 1);

//       // update norm_y
//       double global_dot, local_dot;
//       local_dot = cblas_ddot(dim, y, 1, y, 1);
//       comm.allreduce(1, &local_dot, &global_dot, MPI_SUM);
//       sq_norm_y += global_dot;

//       if (IterativeSolver::os){
// 	*IterativeSolver::os << "Newton: iteration: "
// 			     << iterations << "    "
// 			     << "|p|: " << sqrt(sq_norm_y) << "   "
// 			     << std::endl;
//       }
//     }

//     if(sqrt(sq_norm_y) < tolerance) return true;

//     iterations++;
//   }

//   return false;
// }



// DIRK::LinearOperator implementation
DIRK::LinearOperator::LinearOperator(Communicator &comm, Function &f,
				     const int &dim, 
				     double *&u_tmp, double *&f_tmp) : 
  comm(comm), f(f), dim(dim), u_tmp(u_tmp), f_tmp(f_tmp)
{}


int DIRK::LinearOperator::dim_of_argument(int i) const
{
  return dim;
}


int DIRK::LinearOperator::dim_of_value(int i) const
{
  return dim;
}


void DIRK::LinearOperator::setup(double t, const double *u, double lambda)
{
  this->t = t;
  this->u = u;
  this->lambda = lambda;
}


void DIRK::LinearOperator::operator()(const double *p, double *DFu_p, int i)
{

  double local_dot[2], global_dot[2];
  local_dot[0] = cblas_ddot(dim, u, 1, u, 1);
  local_dot[1] = cblas_ddot(dim, p, 1, p, 1);
  comm.allreduce(2, local_dot, global_dot, MPI_SUM);
  const double norm_u = sqrt(global_dot[0]);
  const double norm_p_sq = global_dot[1];

  const double eps = (norm_p_sq > DBL_EPSILON)?
    sqrt( (1.0+norm_u)*DBL_EPSILON / norm_p_sq ) : sqrt(DBL_EPSILON);
  const double lambda_eps = lambda / eps; 

  dwaxpby(dim, 1.0, u, 1, eps, p, 1, u_tmp, 1);
  f(t, u_tmp, DFu_p);
  for(int i=0; i<dim; i++) DFu_p[i] = p[i] - lambda_eps*(DFu_p[i] - f_tmp[i]);
}





//class DIRK3
static const double delta_dirk = 1.0/2.0 + sqrt(3.0)/6.0;
static const double DIRK3_A[] =
  {delta_dirk, 0.0,
   1.0-2.0*delta_dirk, delta_dirk
  };
static const double DIRK3_b[] = 
  {(0.5-delta_dirk)/(1.0-2.0*delta_dirk), (0.5-delta_dirk)/(1.0-2.0*delta_dirk)};
static const double DIRK3_c[] = 
  {delta_dirk, 1.0-delta_dirk};

DIRK3::DIRK3(Communicator &comm, Function &f) :
  DIRK(comm, 2, 3, f, DIRK3_A, DIRK3_b, DIRK3_c) {}



//class ImplicitEuler
static const double ImplicitEuler_A[] = {1.0};
static const double ImplicitEuler_b[] = {1.0};
static const double ImplicitEuler_c[] = {1.0};

ImplicitEuler::ImplicitEuler(Communicator &comm, Function &f) :
  DIRK(comm, 1, 1, f, ImplicitEuler_A, ImplicitEuler_b, ImplicitEuler_c) 
{}


//class Gauss2
static const double Gauss2_A[] = {0.5};
static const double Gauss2_b[] = {1.0};
static const double Gauss2_c[] = {0.5};

Gauss2::Gauss2(Communicator &comm, Function &f) :
  DIRK(comm, 1, 1, f, Gauss2_A, Gauss2_b, Gauss2_c) 
{}


