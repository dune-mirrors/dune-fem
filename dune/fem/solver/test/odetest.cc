// *****************************
//
// Solve the ODE dy/dt = F(y,t)
//
// *****************************

// standard includes
#include <config.h>
#include <iostream>

// dune includes
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/rungekutta/explicit.hh>


// Faked data structure for our unknown,
// nothing more than a "tuned up" double
class myDest {
  struct SpaceDummy
  {
    int size () const { return 1; }
  };

public:
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;

  myDest(std::string, const SpaceDummy&, const double* u = 0) {
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
  void axpy(RangeFieldType l, const myDest& other) {
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

  RangeFieldType scalarProductDofs ( const myDest &other ) const { return d_ * other.d_; }

private:
  RangeFieldType d_;
};


// implement right hand side F(y,t)
class myRHS : public Dune::Fem::SpaceOperatorInterface<myDest> {
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


int main( int argc, char ** argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // problem data
  const double initialData = 1.0;
  const double startTime = 0.0;
  const double endTime = 2.0;

  // options
  const double stepSize = 0.01;
  const double cfl = 1.;
  const int order = 2;

  // types
  typedef myRHS SpaceOperatorType;
  typedef SpaceOperatorType::DestinationType DestinationType;
  typedef DuneODE::ExplicitRungeKuttaSolver<DestinationType> OdeSolverType;

  // create solver
  typedef Dune::Fem::MPIManager :: CollectiveCommunication  CollectiveCommunication;
  Dune::Fem::TimeProvider<CollectiveCommunication> tp( startTime, cfl, Dune::Fem::MPIManager::comm() );
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

  return 0;
}
