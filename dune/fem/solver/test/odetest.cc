// *****************************
//
// Solve the ODE dy/dt = F(y,t)
//
// *****************************

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/rungekutta.hh>


using namespace Dune;
using namespace DuneODE;
using namespace std;


// Faked data structure for our unknown,
// nothing more than a "tuned up" double
class myDest {
  struct SpaceDummy {};
  
public:
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;
  
  myDest(string, const SpaceDummy&) {
  }
  myDest() {
  }
  RangeFieldType operator[](int i) const {
    return d_;
  }
  RangeFieldType& operator[](int i) {
    return d_;
  }
  void assign(const myDest& other) {
    d_ = other.d_;
  }
  void addScaled(const myDest& other, RangeFieldType l) {
    d_ += other.d_*l;
  }
  myDest& operator*=(RangeFieldType l) {
    d_ *= l;
    return *this;
  }
  myDest& operator+=(const myDest& other) {
    d_ += other.d_;
    return *this;
  }
  myDest& operator-=(const myDest& other) {
    d_ += other.d_;
    return *this;
  }
  
private:
  RangeFieldType d_;
};


// implement right hand side F(y,t)
class myRHS : public SpaceOperatorInterface<myDest> {
public:
  myRHS() {
  }
  
  const SpaceType& space() const {
    return space_;
  }

  void operator()(const DestinationType& x,
                  DestinationType& y) const {
    y[0]=2.0*t_;
    //y[0]=3.0*t_*t_;        
    //y[0]=x[0];
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

  // initialize solution vector
  DestinationType U;
  U[0] = initialData;

  // initialize odesolver
  odeSolver.initialize( U );

  // time loop
  for( tp.init(stepSize); tp.time() < endTime; tp.next(stepSize) ) {
    // do calculation
    odeSolver.solve(U);

    // print out solution
    std::cout << tp.time()
              << " " << U[0]
              << std::endl;
  }
}
