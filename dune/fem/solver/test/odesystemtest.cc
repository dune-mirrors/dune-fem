// ****************************************
//
// Solve the system of ODEs dy/dt = F(y,t)
//
// ****************************************

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/common/fvector.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/rungekutta.hh>


using namespace Dune;
using namespace DuneODE;
using namespace std;


// Data structure for our unknown: fieldvector of
// dimension N (N = number of ODEs) with some additional methods.
// In scalar case (N=1) we have nothing more
// than a simple, tuned up double. 
template <int N>
class myDest : public FieldVector<double, N> {
private:
  struct SpaceDummy {
    int size () const { return N; }
  };
  typedef FieldVector<double, N> BaseType;

public:
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;

  myDest(string, const SpaceDummy&, const double* u = 0) {
  }
  
  myDest() {
  }
  
  void clear() {
    BaseType::operator=(0.);
  }
  
  void assign(const myDest& other) {
    BaseType::operator=(other);
  }
  
  void addScaled(const myDest& other, RangeFieldType l) {
    (*this).axpy(l, other);
  }
  
  double operator()(int i) const {
    if (i<0 || i>this->size-1) {
      std::cout << "ERROR: Accessing element " << i << std::endl;
    }
    return (*this)[i];
  }
};


// implement right hand side F(y,t)
// here: system of three ODEs
class myRHS : public SpaceOperatorInterface<myDest<3> > {
public:
  myRHS() {
  }
  
  const SpaceType& space() const {
    return space_;
  }

  void operator()(const DestinationType& x,
                  DestinationType& y) const {
    y[0]=2.0*t_;
    y[1]=3.0*t_*t_;
    y[2]=x[2];
  }

  void setTime(const double time) {
    t_=time;
  }

private:
  SpaceType space_;
  double t_;
};


int main() {
  // problem data
  const double initialData = 1.0;
  const double startTime = -2.0;
  const double endTime = 2.0;

  // options
  const double stepSize = 0.01;
  const double cfl = 1.;
  const int order = 2;

  // types
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;
  typedef ExplicitRungeKuttaSolver<DestinationType> OdeSolverType;

  // create solver
  TimeProvider<> tp( startTime, cfl );
  SpaceOperatorType spaceOperator;
  OdeSolverType odeSolver( spaceOperator, tp, order );

  // initialize solution vector, same initial data for all components
  DestinationType U;
  for (int i=0; i<U.size; i++)
    U[i] = initialData;

  // initialize odesolver
  odeSolver.initialize( U );

  // time loop
  for( tp.init(stepSize); tp.time() < endTime; tp.next(stepSize) ) {
    // do calculation
    odeSolver.solve(U);

    // print out solution
    std::cout << tp.time()
              << " " << U
              << std::endl;
  }
}
