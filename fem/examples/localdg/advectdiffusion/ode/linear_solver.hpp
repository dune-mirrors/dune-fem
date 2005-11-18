#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include <iostream>
#include <cmath>
#include "communicator.hpp"
#include "dynamical_object.hpp"
#include "function.hpp"
#include "iterative_solver.hpp"
#include "matrix.hpp"


namespace DuneODE {

class DirectLinearSolver
{
public:
  DirectLinearSolver();

  // decomposition of matrix a, a can be modified
  virtual void prepare(int dim, double *a);

  // solve A x = b after A is "prepared", b is initially stored in x
  virtual bool solve(double *x) = 0;

protected:
  int dim;
  double *a; // pointer to matrix
};



class IterativeLinearSolver : public IterativeSolver
{
public:
  IterativeLinearSolver(const Communicator &comm);
  virtual ~IterativeLinearSolver();

  virtual void set_preconditioner(Function &preconditioner);

  // solve Au = b,   A = linear_operator
  // return convergence
  virtual bool solve(Function &op, double *u, const double *b) = 0;

protected:
  // Communicator
  const Communicator &comm;

  int dim;

  // we want to solve Au = b, the preconditioner gives z = M^{-1} v,
  // where M is an approximation to A
  Function *preconditioner;
};





class LUSolver : public DirectLinearSolver, public DynamicalObject
{
public:
  LUSolver();
  ~LUSolver();

  // from DirectLinearSolver
  // build the LU dcomposition
  virtual void prepare(int dim, double *a);
  virtual bool solve(double *x);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  int *p; // permutation array
};


class QRSolver : public DirectLinearSolver, public DynamicalObject
{
public:
  QRSolver();
  ~QRSolver();

  // from DirectLinearSolver
  // build the QR dcomposition
  virtual void prepare(int dim, double *a);
  virtual bool solve(double *x);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  double *d; // diagonal elements of the R-matrix
};




// Saad, Youcef;  Schultz, Martin H.
// GMRES: A generalized minimal residual algorithm for solving nonsymmetric 
// linear systems. (English)
// [J] SIAM J. Sci. Stat. Comput. 7, 856-869 (1986). [ISSN 0196-5204]
template<int m>
class GMRES : public IterativeLinearSolver, public DynamicalObject
{
public:
  GMRES(const Communicator &comm);
  virtual ~GMRES();
  virtual void set_preconditioner(Function &preconditioner);

  // from Function, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  Matrix H; // \in \R^{m+1 \times m}
  double g[m+1], s[m], c[m], y[m], local_dot[m], global_dot[m], *v[m+1], *z;
};



// Saad, Youcef
// A flexible inner-outer preconditioned GMRES algorithm. (English)
// [J] SIAM J. Sci. Comput. 14, No.2, 461-469 (1993). [ISSN 1064-8275]
template<int m>
class FGMRES : public IterativeLinearSolver, public DynamicalObject
{
public:
  FGMRES(const Communicator &comm);
  virtual ~FGMRES();

  // solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  Matrix H; // \in \R^{m+1 \times m}
  double g[m+1], s[m], c[m], y[m], local_dot[m], global_dot[m], *v[m+1], *z[m];
};




// Iterative methods for sparse linear systems / Yousef Saad. 
// Boston, Mass. : PWS Publ., 1996. - XVI, 447 S. : graph. Darst.; 
// (engl.)(The PWS series in computer science)
// ISBN 0-534-94776-X 
//
// BICGSTAB - Solver page 220
class BICGSTAB : public IterativeLinearSolver, public DynamicalObject
{
public:
  BICGSTAB(const Communicator &comm);
  ~BICGSTAB();

  // from IterativeLinearSolver, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  double *r, *r_star, *p, *s, *tmp;
};




// CG and Preconditioned CG scheme
// Braess - Finite Elemente
//
// CG/PCG - Solver
class CG : public IterativeLinearSolver, public DynamicalObject
{
public:
  CG(const Communicator &comm);
  virtual ~CG();

  // from IterativeLinearSolver, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);


private:
  double *r, *d, *h;
};





// inline implementations

// IterativeLinearSolver
inline
IterativeLinearSolver::IterativeLinearSolver(const Communicator &comm) :
  comm(comm), dim(0), preconditioner(NULL) {}

inline
IterativeLinearSolver::~IterativeLinearSolver() {}

inline
void IterativeLinearSolver::set_preconditioner(Function &preconditioner)
{
  this->preconditioner = &preconditioner;
}



// DirectLinearSolver
inline DirectLinearSolver::DirectLinearSolver() : dim(0), a(NULL) {}

inline
void DirectLinearSolver::prepare(int dim, double *a)
{
  this->dim = dim;
  this->a = a;
}


// template solvers
// GMRES inline implementation

#include <cmath>
#include <cassert>
#include "blas.hpp"


template<int m>
inline
GMRES<m>::GMRES(const Communicator &comm) : 
  IterativeLinearSolver(comm), DynamicalObject("GMRES", comm.id()), 
  H(m+1,m), z(NULL)
{
  for(int i=0; i<=m; i++) v[i] = NULL;
}



template<int m>
inline
GMRES<m>::~GMRES()
{
  delete[] v[0];
  delete[] z;
}



template<int m>
inline
void GMRES<m>::set_preconditioner(Function &preconditioner)
{
  this->preconditioner = &preconditioner;
  if (dim > 0){
    z = new double[size()];
    assert(z);
  }
}


template<int m>
inline
void GMRES<m>::resize(int new_size, int component)
{
  if (preconditioner){
    delete[] z;
    z = new double[new_size];
    assert(z);
  }

  delete[] v[0];
  v[0] = new double[(m+1)*new_size];
  assert(v[0]);
  for(int i=1; i<=m; i++) v[i] = v[0] + i*new_size;
}
 



template<int m>
inline
bool GMRES<m>::solve(Function &op, double *u, const double *b)
{
  dim = op.dim_of_value();
  new_size(dim);

  int iterations = 0;
  while (true){
    // start
    op(u, v[0]);
    cblas_daxpy(dim, -1.0, b, 1, v[0], 1);
    local_dot[0] = cblas_ddot(dim, v[0], 1, v[0], 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    double res = sqrt(global_dot[0]); 
    if (res < tolerance) break;
    g[0] = -res;
    for(int i=1; i<=m; i++) g[i] = 0.0;
    cblas_dscal(dim, 1.0/res, v[0], 1);
    
    // iterate
    for(int j=0; j<m; j++){
      // apply the linear operator (perhaps in combination with the 
      // preconditioner)
      if (preconditioner){
	(*preconditioner)(v[j], z);
	op(z, v[j+1]);
      }
      else op(v[j], v[j+1]);

      for(int i=0; i<=j; i++) local_dot[i]=cblas_ddot(dim, v[j+1], 1, v[i], 1);
      comm.allreduce(j+1, local_dot, global_dot, MPI_SUM);      
      for(int i=0; i<=j; i++) H(i,j) = global_dot[i]; 

      for(int i=0; i<=j; i++) cblas_daxpy(dim, -H(i,j), v[i], 1, v[j+1], 1);
      local_dot[0] = cblas_ddot(dim, v[j+1], 1, v[j+1], 1);
      comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
      H(j+1,j) = sqrt(global_dot[0]); 
      cblas_dscal(dim, 1.0/H(j+1,j), v[j+1], 1);

      // perform Givens rotation
      for(int i=0; i<j; i++){
	cblas_drot(1, &H(i+1,j), 1, &H(i,j), 1, c[i], s[i]);
      }
      const double h_j_j = H(j,j);
      const double h_jp_j = H(j+1,j);
      const double norm = sqrt(h_j_j*h_j_j + h_jp_j*h_jp_j);
      c[j] = h_j_j / norm;
      s[j] = -h_jp_j / norm;
      cblas_drot(1, &H(j+1,j), 1, &H(j,j), 1, c[j], s[j]);
      cblas_drot(1, &g[j+1], 1, &g[j], 1, c[j], s[j]);

      //*os << fabs(g[j+1]) << std::endl;

      iterations++;
      if (fabs(g[j+1]) < tolerance 
	  || iterations >= max_num_of_iterations) break;
    }
    
    //
    // form the approximate solution
    //

    int last = iterations%m;
    if (last == 0) last = m;
    // compute y via backsubstitution
    for(int i=last-1; i>=0; i--){
      const double dot = cblas_ddot(last-(i+1), &H(i,i)+1, 1, &y[i+1], 1);
      y[i] = (g[i] - dot)/ H(i,i);
    }    

    // update the approx. solution
    if (preconditioner){
      // u += M^{-1} (v[0], ..., v[last-1]) y	
      double *u_tmp = v[m]; // we don't need this vector anymore
      dset(dim, 0.0, u_tmp, 1);
      for(int i=0; i<last; i++)	cblas_daxpy(dim, y[i], v[i], 1, u_tmp, 1);
      (*preconditioner)(u_tmp, z);
      cblas_daxpy(dim, 1.0, z, 1, u, 1);
    }
    else{
      // u += (v[0], ..., v[last-1]) y
      for(int i=0; i<last; i++)	cblas_daxpy(dim, y[i], v[i], 1, u, 1);
    }

    if (fabs(g[last]) < tolerance) break;
  }

  // output
  if (IterativeSolver::os){
    *IterativeSolver::os << "GMRES " << comm.id() 
			 << ": number of iterations: "      
			 << iterations
			 << std::endl;
  }

  return (iterations < max_num_of_iterations)? true: false;
}

// FGMRES inline implementation

#include <cmath>
#include <cassert>
#include "blas.hpp"


template<int m>
inline
FGMRES<m>::FGMRES(const Communicator &comm) : 
  IterativeLinearSolver(comm), DynamicalObject("FGMRES", comm.id()), H(m+1,m)
{
  for(int i=0; i<m; i++) v[i] = z[i] = NULL;
  v[m] = NULL;
}


template<int m>
inline
FGMRES<m>::~FGMRES()
{
  delete[] v[0];
}


template<int m>
inline
void FGMRES<m>::resize(int new_size, int component)
{
  delete[] v[0];
  v[0] = new double[(2*m+1)*new_size];
  assert(v[0]);
  for(int i=0; i<m; i++){
    v[i] = v[0] + i*new_size;
    z[i] = v[0] + (m+1 + i)*new_size;
  }
  v[m] = v[0] + m*new_size;
}


template<int m>
inline
bool FGMRES<m>::solve(Function &op, double *u, const double *b)
{
  // it does not make sense to use FGMRES without a preconditioner
  assert(preconditioner);
  dim = op.dim_of_value();
  new_size(dim);

  int iterations = 0;
  while (true){
    // start
    op(u, v[0]);
    cblas_daxpy(dim, -1.0, b, 1, v[0], 1);
    local_dot[0] = cblas_ddot(dim, v[0], 1, v[0], 1);
    comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
    double res = sqrt(global_dot[0]); 
    if (res < tolerance) break;
    g[0] = -res;
    for(int i=1; i<=m; i++) g[i] = 0.0;
    cblas_dscal(dim, 1.0/res, v[0], 1);

    // iterate
    for(int j=0; j<m; j++){
      // apply preconditioner and the linear operator
      (*preconditioner)(v[j], z[j]);
      op(z[j], v[j+1]);      
      
      // Gram-Schmidt procedure
      for(int i=0; i<=j; i++) local_dot[i]=cblas_ddot(dim, v[j+1], 1, v[i], 1);
      comm.allreduce(j+1, local_dot, global_dot, MPI_SUM);      
      for(int i=0; i<=j; i++) H(i,j) = global_dot[i]; 

      for(int i=0; i<=j; i++) cblas_daxpy(dim, -H(i,j), v[i], 1, v[j+1], 1);
      local_dot[0] = cblas_ddot(dim, v[j+1], 1, v[j+1], 1);
      comm.allreduce(1, local_dot, global_dot, MPI_SUM);      
      H(j+1,j) = sqrt(global_dot[0]); 
      cblas_dscal(dim, 1.0/H(j+1,j), v[j+1], 1);

      // perform Givens rotation
      for(int i=0; i<j; i++){
        cblas_drot(1, &H(i+1,j), 1, &H(i,j), 1, c[i], s[i]);
      }
      const double h_j_j = H(j,j);
      const double h_jp_j = H(j+1,j);
      const double norm = sqrt(h_j_j*h_j_j + h_jp_j*h_jp_j);
      c[j] = h_j_j / norm;
      s[j] = -h_jp_j / norm;
      cblas_drot(1, &H(j+1,j), 1, &H(j,j), 1, c[j], s[j]);
      cblas_drot(1, &g[j+1], 1, &g[j], 1, c[j], s[j]);

      //*output_stream << fabs(g[j+1]) << std::endl;

      iterations++;
      if (fabs(g[j+1]) < tolerance 
	  || iterations >= max_num_of_iterations) break;
    }

    //
    // form the approximate solution
    //
      
    int last = iterations%m;
    if (last == 0) last = m;
    // compute y via backsubstitution
    for(int i=last-1; i>=0; i--){
      const double dot = cblas_ddot(last-(i+1), &H(i,i)+1, 1, &y[i+1], 1);
      y[i] = (g[i] - dot) / H(i,i);
    }    

    // update the approx. solution
    // u += (z[0], ..., z[last-1]) y
    for(int i=0; i<last; i++) cblas_daxpy(dim, y[i], z[i], 1, u, 1);
    
    if (fabs(g[last]) < tolerance) break;
  }

  // output
  if (IterativeSolver::os){
    *IterativeSolver::os<< "FGMRES " << comm.id() 
			<< ": number of iterations: "
			<< iterations
			<< std::endl;
  }

  return (iterations < max_num_of_iterations)? true: false;
}

 
}



#endif

