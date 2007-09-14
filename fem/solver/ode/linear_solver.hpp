#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include <iostream>
#include <cmath>
#include "communicator.hpp"
#include "dynamical_object.hpp"
#include "function.hpp"
#include "iterative_solver.hpp"
#include "matrix.hpp"

namespace pardg
{


class DirectLinearSolver
{
public:
  DirectLinearSolver();
  virtual ~DirectLinearSolver() {}

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
  IterativeLinearSolver(Communicator &comm);
  virtual ~IterativeLinearSolver();

  virtual void set_preconditioner(Function &preconditioner);

  // set pointer to predonditioner to zero 
  virtual void unset_preconditioner() 
  {
    preconditioner = 0;
  }
  
  // solve Au = b,   A = linear_operator
  // return convergence
  virtual bool solve(Function &op, double *u, const double *b) = 0;

protected:
  // Communicator
  Communicator &comm;

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
class GMRES : public IterativeLinearSolver, public DynamicalObject
{
public:
  GMRES(Communicator &comm, int m);
  virtual ~GMRES();
  virtual void set_preconditioner(Function &preconditioner);
  // set pointer to predonditioner to zero 
  virtual void unset_preconditioner();

  // from Function, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  const int m;
  Matrix H; // \in \R^{m+1 \times m}
  double *g, *s, *c, *y, *local_dot, *global_dot, *v, *z;
};



// Saad, Youcef
// A flexible inner-outer preconditioned GMRES algorithm. (English)
// [J] SIAM J. Sci. Comput. 14, No.2, 461-469 (1993). [ISSN 1064-8275]
class FGMRES : public IterativeLinearSolver, public DynamicalObject
{
public:
  FGMRES(Communicator &comm, int m);
  virtual ~FGMRES();

  // solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  const int m;
  Matrix H; // \in \R^{m+1 \times m}
  double *g, *s, *c, *y, *local_dot, *global_dot, *v, *z;
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
  BICGSTAB(Communicator &comm);
  ~BICGSTAB();

  virtual void set_preconditioner(Function &preconditioner);
  virtual void unset_preconditioner();

  // from IterativeLinearSolver, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);
  virtual bool solve_old(Function &linear_operator,double *u,const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  double *r, *r_star, *p, *s, *tmp, *z;
};




// CG and Preconditioned CG scheme
// Braess - Finite Elemente
//
// CG/PCG - Solver
class CG : public IterativeLinearSolver, public DynamicalObject
{
public:
  CG(Communicator &comm);
  virtual ~CG();

  // from IterativeLinearSolver, solve Au = b, Au = linear_operator(u)
  virtual bool solve(Function &linear_operator, double *u, const double *b);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);


private:
  double *r, *d, *h;
};



} // namespace pardg



// inline implementations

// IterativeLinearSolver
inline
pardg::IterativeLinearSolver::IterativeLinearSolver(Communicator &comm) :
  comm(comm), dim(0), preconditioner(NULL) {}

inline
pardg::IterativeLinearSolver::~IterativeLinearSolver() {}

inline
void pardg::IterativeLinearSolver::set_preconditioner(Function &preconditioner)
{
  this->preconditioner = &preconditioner;
}



// DirectLinearSolver
inline pardg::DirectLinearSolver::DirectLinearSolver() : dim(0), a(NULL) {}

inline
void pardg::DirectLinearSolver::prepare(int dim, double *a)
{
  this->dim = dim;
  this->a = a;
}




#endif

