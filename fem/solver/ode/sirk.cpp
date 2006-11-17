#include <cfloat>
#include "ode_solver.hpp"

using namespace pardg;



SIRK::SIRK(Communicator &comm, int num_of_stages, int order, 
	   Function &f, Function &fex,
	   const double *a, const double *b, const double *c,
	   const double *aex, const double *cex) :
  ODESolver(comm, 0), 
  f(f), fex(fex), 
  num_of_stages(num_of_stages), order(order),
  A(num_of_stages, a), b(num_of_stages, b), c(num_of_stages, c),
  Aex(num_of_stages, aex), cex(num_of_stages, cex),
  alpha(num_of_stages), beta(num_of_stages), gamma(num_of_stages),
  alphaex(num_of_stages), gammaex(num_of_stages),
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

  // build matrix alphaex
  alphaex = Aex * Ainv;
  
  // build vector beta
  for(int i=0; i<num_of_stages; i++){
    beta[i] = 0.0;
    for(int j=0; j<num_of_stages; j++) beta[i] += b[j] * Ainv(j,i);
  }
  
  // build vector gamma / gammaex
  for(int i=0; i<num_of_stages; i++){
    gamma[i] = 1.0;
    gammaex[i] = 1.0;
    for(int j=0; j<i; j++){
      gamma[i] -= alpha(i,j);
      gammaex[i] -= alphaex(i,j);
    }
  }

  // delta
  delta = 1.0;
  for(int i=0; i<num_of_stages; i++) delta -= beta[i];
}



void SIRK::set_linear_solver(IterativeLinearSolver &ls)
{
  ils = &ls;
}


void SIRK::resize(int new_size, int component)
{
  // new_size >= dim
  delete[] U;
  U = new double[ (num_of_stages + 5) * new_size ];  
  Fpre = U + num_of_stages * new_size;
  F = Fpre + new_size;    
  y = F + new_size;
  u_tmp = y + new_size;    // for LinearOperator
  f_tmp = u_tmp + new_size; // 
}



bool SIRK::step(double t, double dt, double *u)
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


bool SIRK::step_iterative(double t, double dt, double *u)
{
  for(int i=0; i<num_of_stages; i++){
    double *ui = U+i*dim;

    // setup of Fpre, store ui_ex in ui and f(ui_ex) in y
    const double _gamma = gamma[i];
    const double _gammaex = gammaex[i];
    for(int l=0; l<dim; l++){
      Fpre[l] = _gamma * u[l];
      ui[l] = _gammaex * u[l];
    }
    for(int j=0; j<i; j++){
      cblas_daxpy(dim, alpha(i,j), U+j*dim, 1, Fpre, 1);
      cblas_daxpy(dim, alphaex(i,j), U+j*dim, 1, ui, 1);
    }
    fex(t + cex[i]*dt, ui, y);
    cblas_daxpy(dim, alpha(i,i)*dt, y, 1, Fpre, 1);

    // prediction uf ui, todo: extrapolation or something...
    // ui = u^n
    cblas_dcopy(dim, u, 1, ui, 1);

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


// bool SIRK::step_iterative(double t, double dt, double *u)
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



// SIRK::LinearOperator implementation
SIRK::LinearOperator::LinearOperator(Communicator &comm, Function &f,
				     const int &dim, 
				     double *&u_tmp, double *&f_tmp) : 
  comm(comm), f(f), dim(dim), u_tmp(u_tmp), f_tmp(f_tmp)
{}


int SIRK::LinearOperator::dim_of_argument(int i) const
{
  return dim;
}


int SIRK::LinearOperator::dim_of_value(int i) const
{
  return dim;
}


void SIRK::LinearOperator::setup(double t, const double *u, double lambda)
{
  this->t = t;
  this->u = u;
  this->lambda = lambda;
}


void SIRK::LinearOperator::operator()(const double *p, double *DFu_p, int i)
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





//class SemiImplicitEuler, 1 stage, 1st order
// implicit part
static const double SIEuler_A[] = {1.0};
static const double SIEuler_c[] = {1.0};
// explicit part
static const double SIEuler_Aex[] = {0.0};
static const double SIEuler_cex[] = {0.0};
// common part
static const double SIEuler_b[] = {1.0};

SemiImplicitEuler::SemiImplicitEuler(Communicator &comm, 
				     Function &f, Function &fex) :
  SIRK(comm, 1, 1, f, fex, SIEuler_A, SIEuler_b, SIEuler_c, 
       SIEuler_Aex, SIEuler_cex) 
{}





//class SIRK23, 3 stages, 2nd order
// implicit part
static const double SIRK23_A[] =
  {0.5, 0.0, 0.0,
   -1.0, 0.5, 0.0,
   0.25, 0.25, 0.5
  };
static const double SIRK23_c[] = 
  {0.5, -0.5, 1.0};
// explicit part
static const double SIRK23_Aex[] =
  {0.0, 0.0, 0.0,
   1.0, 0.0, 0.0,
   0.0, 0.5, 0.0
  };
static const double SIRK23_cex[] = 
  {0.0, 1.0, 0.5};
// common part
static const double SIRK23_b[] = 
  {0.25, 0.25, 0.5};

SIRK23::SIRK23(Communicator &comm, Function &f, Function &fex) :
  SIRK(comm, 3, 2, f, fex, SIRK23_A, SIRK23_b, SIRK23_c, 
       SIRK23_Aex, SIRK23_cex) 
{}



//class SIRK33, 3 stages, 3rd order
// implicit part
static const double SIRK33_A[] =
  {3.0/4.0, 0.0, 0.0,
   5589.0/6524.0, 75.0/233.0, 0.0,
   7691.0/26096.0, -26335.0/78288.0, 65.0/168.0
  };
static const double SIRK33_c[] = 
  {3.0/4.0, 
   5589.0/6524.0 + 75.0/233.0, 
   7691.0/26096.0 - 26335.0/78288.0 + 65.0/168.0
  };
// explicit part
static const double SIRK33_Aex[] =
  {0.0, 0.0, 0.0,
   8.0/7.0, 0.0, 0.0,
   71.0/252.0, 7.0/36.0, 0.0
  };
static const double SIRK33_cex[] = 
  {0.0, 8.0/7.0, 71.0/252.0 + 7.0/36.0};
// common part
static const double SIRK33_b[] = 
  {1.0/8.0, 1.0/8.0, 3.0/4.0};

SIRK33::SIRK33(Communicator &comm, Function &f, Function &fex) :
  SIRK(comm, 3, 2, f, fex, SIRK33_A, SIRK33_b, SIRK33_c, 
       SIRK33_Aex, SIRK33_cex) 
{}



//class IMEX_SSP222, 2 stages, 2nd order
// implicit part
static const double delta = 1.0 - 1.0/sqrt(2.0);
static const double IMEX_SSP222_A[] =
  {delta, 0.0,
   1.0-2.0*delta, delta
  };
static const double IMEX_SSP222_c[] = 
  {delta, 1.0-delta};
// explicit part
static const double IMEX_SSP222_Aex[] =
  {0.0, 0.0,
   1.0, 0.0
  };
static const double IMEX_SSP222_cex[] = 
  {0.0, 1.0};
// common part
static const double IMEX_SSP222_b[] = 
  {0.5, 0.5};

IMEX_SSP222::IMEX_SSP222(Communicator &comm, 
			 Function &f, Function &fex) :
  SIRK(comm, 2, 2, f, fex, IMEX_SSP222_A, IMEX_SSP222_b, IMEX_SSP222_c, 
       IMEX_SSP222_Aex, IMEX_SSP222_cex) 
{}



