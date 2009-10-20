#include <cmath> 
#include <cassert>
#include <cfloat>
#include "newton.hpp"
#include "blas.hpp"
#include "matrix_decomposition.hpp"



Newton::Newton() : 
  DynamicalObject("Newton", -1), f(NULL), Df(NULL), p(NULL)
{
  // set this to some useful values
  tolerance = 1.0e-6;
  max_num_of_iterations = 20;
}


Newton::~Newton()
{
  delete[] f;
  delete[] p;
}


void Newton::resize(int new_size, int component)
{
  delete[] f;
  delete[] p;
  f = new double[(1 + new_size) * new_size]; // new_size >= dim
  p = new int[new_size];
  assert(f && p);
  Df = f + new_size;
}
 

bool Newton::solve(Function &F, double *u)
{
  dim = F.dim_of_argument();
  new_size(dim);
  int num_of_iterations = 0;
  double res_old = DBL_MAX;

  while (num_of_iterations < max_num_of_iterations) {
    F(u, f);
    
    const double res = cblas_dnrm2(dim, f, 1); 
    if (res < tolerance) return true;
    if (res >= res_old) return false;
    res_old = res;

    F(u, Df, 1);
    LU_decomposition(dim, Df, p);
    LU_solve(dim, Df, p, f);
        
    cblas_daxpy(dim, -1.0, f, 1, u, 1);

    if (NonlinearSolver::os){    
      *NonlinearSolver::os << "Newton: iteration: " 
			   << num_of_iterations << "    "
			   << "res: " << res << "   "
			   << std::endl;
    }

    num_of_iterations++;
  } 
  
  return false; // num_of_iterations >= max_num_of_iterations
}


