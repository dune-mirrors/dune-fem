#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <iostream>
#include "communicator.hpp"
#include "function.hpp"
#include "limiter.hpp"
#include "dynamical_object.hpp"
#include "iterative_solver.hpp"
#include "linear_solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"


namespace pardg
{

class ODESolver : public DynamicalObject
{
public:
  ODESolver(Communicator &comm, int num_of_tmpobj);
  virtual ~ODESolver();

  void set_output(std::ostream &os);
  void set_limiter(Limiter &limiter);

  // user-interface for solving
  virtual bool step(double t, double dt, double *u) = 0; 

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

  Communicator &comm;
  const int num_of_tmpobj;
  int dim;
  double *U;
  Limiter *limiter;
  std::ostream *os;
};



// Diagonally implicit Runge Kutta schemes
class DIRK : public ODESolver, public IterativeSolver
{
public:
  DIRK(Communicator &comm, int num_of_stages, int order, Function &f,
       const double *a, const double *b, const double *c);

  void set_linear_solver(IterativeLinearSolver &ls);

  // from ODESolver
  virtual bool step(double t, double dt, double *u);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);
  
private:
  class LinearOperator : public Function
  {
  public:
    LinearOperator(Communicator &comm, Function &f, const int &dim,
		   double *&u_tmp, double *&f_tmp);

    void setup(double t, const double *u, double lambda);

    // from Function
    virtual void operator()(const double *p, double *DFu_p, int i = 0);
    virtual int dim_of_argument(int i = 0) const;
    virtual int dim_of_value(int i = 0) const;

  private:
    Communicator &comm;
    Function &f;
    const int &dim;
    double *&u_tmp, *&f_tmp;
    double lambda, t;
    const double *u;
  };


  Function &f;
  const int num_of_stages;
  const int order;
  const Matrix A;    // classical Butcher table
  const Vector b,c;  //
  Matrix alpha; // modified table
  Vector beta, gamma;
  double delta;
  double *F, *Fpre, *y;

  // iterative solver
  bool step_iterative(double t, double dt, double *u);
  IterativeLinearSolver *ils;
  LinearOperator op;
  double *u_tmp, *f_tmp; 
};



// semi implicit Runge Kutta schemes
class SIRK : public ODESolver, public IterativeSolver
{
public:
  SIRK(Communicator &comm, int num_of_stages, int order, 
       Function &f, Function &fex,
       const double *a, const double *b, const double *c,
       const double *aex, const double *cex);

  void set_linear_solver(IterativeLinearSolver &ls);

  // from ODESolver
  virtual bool step(double t, double dt, double *u);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);
  
private:
  class LinearOperator : public Function
  {
  public:
    LinearOperator(Communicator &comm, Function &f, const int &dim,
		   double *&u_tmp, double *&f_tmp);

    void setup(double t, const double *u, double lambda);

    // from Function
    virtual void operator()(const double *p, double *DFu_p, int i = 0);
    virtual int dim_of_argument(int i = 0) const;
    virtual int dim_of_value(int i = 0) const;

  private:
    Communicator &comm;
    Function &f;
    const int &dim;
    double *&u_tmp, *&f_tmp;
    double lambda, t;
    const double *u;
  };


  Function &f, &fex; // implicit / explicit part
  const int num_of_stages;
  const int order;
  const Matrix A;    // classical Butcher table
  const Vector b,c;  // for implicit part
  const Matrix Aex;  // classical Butcher table
  const Vector cex;  // for explicit part
  Matrix alpha;      // modified table for implicit part
  Vector beta, gamma;//
  Matrix alphaex;    // modified table for explicit part
  Vector gammaex;    //
  double delta;
  double *F, *Fpre, *y;

  // iterative solver
  bool step_iterative(double t, double dt, double *u);
  IterativeLinearSolver *ils;
  LinearOperator op;
  double *u_tmp, *f_tmp; 
};




class ExplicitSSP : public ODESolver
{
public:
  ExplicitSSP(Communicator &comm, Function &f, int num_of_stages);

  // from ODESolver
  virtual bool step(double t, double dt, double *u);
  
private:
  Function &f;
  const int num_of_stages;
  const Vector alpha;  
};




class ExplicitBulirschStoer : public ODESolver
{
public:
  ExplicitBulirschStoer(Communicator &comm, Function &f, 
			int num_of_stages, 
			int (*seq)(int) = DoubleHarmonicSequence );

  virtual bool step(double t, double dt, double *u);

  static int DoubleRombergSequence(int i);
  static int DoubleHarmonicSequence(int i);

private:
  const int num_of_stages;
  Function &f;
  int (*sequence)(int i);
};




class ImplicitBulirschStoer : public ODESolver, public IterativeSolver
{
public:
  ImplicitBulirschStoer(Communicator &comm, Function &f, 
			int num_of_stages, 
			int (*seq)(int) = DoubleHarmonicSequence );

  virtual bool step(double t, double dt, double *u);
  void set_linear_solver(DirectLinearSolver &ls);
  void set_linear_solver(IterativeLinearSolver &ls);

  static int DoubleRombergSequence(int i);
  static int DoubleHarmonicSequence(int i);

protected:
  // from DynamicalObject
  virtual void resize(int new_size, int component);

private:
  class LinearOperator : public Function
  {
  public:
    LinearOperator(ImplicitBulirschStoer &ibs);

    // from Function
    virtual void operator()(const double *p, double *DFu_p, int i = 0);
    virtual int dim_of_argument(int i = 0) const;
    virtual int dim_of_value(int i = 0) const;

  private:
    ImplicitBulirschStoer &ibs;
  };

  // common stuff
  const int num_of_stages;
  Function &f;
  int (*sequence)(int i);
  double *z_k, *z_km, *F, *tmp;
  double t, dt_n;
  int k;

  // direct solver
  bool step_direct(double t, double dt, double *u);
  DirectLinearSolver *dls;

  // iterative solver
  bool step_iterative(double t, double dt, double *u);
  IterativeLinearSolver *ils;
  LinearOperator op;
};






class ExplicitRungeKutta : public ODESolver
{
public:
  virtual bool step(double t, double dt, double *u);
  
protected:
  ExplicitRungeKutta(Communicator &comm, 
		     int num_of_stages, int order, Function &f);
  ExplicitRungeKutta(Communicator &comm, 
		     int num_of_stages, int order, Function &f,
		     const double *a, const double *b, const double *c);

  // solver interface 
  virtual bool solve(double t, double dt, double *u); 

  Function &f;
  const int num_of_stages;

private:
  const Matrix A;    // classical Butcher table
  const Vector b,c;  //
  Matrix alpha;      // "modified Butcher" table
  Vector beta, gamma;//
  const int order;
};



class ExplicitEuler : public ExplicitRungeKutta
{
public:
  ExplicitEuler(Communicator &comm, Function &f);
  virtual bool solve(double t, double dt, double *u);
};



class ExplicitModifiedEuler : public ExplicitRungeKutta
{
public:
  ExplicitModifiedEuler(Communicator &comm, Function &f);
  virtual bool solve(double t, double dt, double *u);
};



class ExplicitTVD2 : public ExplicitRungeKutta
{
public:
  ExplicitTVD2(Communicator &comm, Function &f);
  virtual bool solve(double t, double dt, double *u);
};




class ExplicitTVD3 : public ExplicitRungeKutta
{
public:
  ExplicitTVD3(Communicator &comm, Function &f);
  virtual bool solve(double t, double dt, double *u);
};



class ExplicitRK3 : public ExplicitRungeKutta
{
public:
  ExplicitRK3(Communicator &comm, Function &f);
  //virtual bool solve(double t, double dt, double *u);
};



class ExplicitRK4 : public ExplicitRungeKutta
{
public:
  ExplicitRK4(Communicator &comm, Function &f);
  virtual bool solve(double t, double dt, double *u);
};


class ExplicitButcher6 : public ExplicitRungeKutta
{
public:
  ExplicitButcher6(Communicator &comm, Function &f);
  //virtual bool solve(double t, double dt, double *u);
};



class ExplicitRK4b : public ExplicitRungeKutta
{
public:
  ExplicitRK4b(Communicator &comm, Function &f);
};



// DIRK methods
class ImplicitEuler : public DIRK
{
public:
  ImplicitEuler(Communicator &comm, Function &f);
};

class Gauss2 : public DIRK
{
public:
  Gauss2(Communicator &comm, Function &f);
};

class DIRK3 : public DIRK
{
public:
  DIRK3(Communicator &comm, Function &f);
};




// SIRK methods
class SemiImplicitEuler : public SIRK
{
public:
  SemiImplicitEuler(Communicator &comm, Function &f, Function &fex);
};


class SIRK23 : public SIRK
{
public:
  SIRK23(Communicator &comm, Function &f, Function &fex);
};


class SIRK33 : public SIRK
{
public:
  SIRK33(Communicator &comm, Function &f, Function &fex);
};


class IMEX_SSP222 : public SIRK
{
public:
  IMEX_SSP222(Communicator &comm, Function &f, Function &fex);
};



} // namespace pardg



#endif
